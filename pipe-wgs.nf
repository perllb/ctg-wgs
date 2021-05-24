#!/usr/bin/env nextFlow

// set variables
runfolder = params.runfolder
basedir = params.basedir
metaID = params.metaid
OUTDIR = params.outdir
FQDIR = params.fqdir
QCDIR = params.qcdir
CTGQC = params.ctgqc
demux = params.demux
b2farg = params.bcl2fastq_arg
index = params.index

// Read and process sample sheet
sheet = file(params.sheet)

// samplesheet to be parsed as input channel (take everything below [Data])
channel_sheet = file("$basedir/samplesheet.channel.nf.ctg-wgs.csv")

// create new samplesheet parsed to fit the format for dragen demux
newsheet = "${basedir}/samplesheet.demux.nf.ctg-wgs.csv"

// Read and process sample sheet
all_lines = sheet.readLines()
write_b = false // if next lines has sample info
channel_sheet.text=""     

for ( line in all_lines ) {

    if ( write_b ) {
	channel_sheet.append(line + "\n")
    }
    if (line.contains("[Data]")) {
	write_b = true
    }
}


println "============================="
println ">>> exome dragen pipeline "
println ""
println "> INPUT: "
println "> runfolder		: $runfolder "
println "> sample-sheet		: $sheet "
println "> run-meta-id		: $metaID "
println "> basedir		: $basedir "
println "> bcl2fastq args	: $b2farg "
println ""
println "> OUTPUT: "
println "> output-dir		: $OUTDIR "
println "> fastq-dir		: $FQDIR "
println "> qc-dir		: $QCDIR "
println "> ctg-qc-dir		: $CTGQC "
println "============================="

// sample info
Channel
    .fromPath(channel_sheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.Sample_ID, row.Sample_Project, row.Sample_Species) }
    .unique()
    .tap{infoSamples}
    .set{ move_fastq_csv; }

println " > Samples to process: "
infoSamples.subscribe{ println "Sample: $it" }

// Parse samplesheet
process parsesheet {

	tag "$metaID"

	input:
	val sheet
	val index

	output:
	val newsheet into demux_sheet

	when:
	demux == 'y'

	"""
#!/opt/conda/bin/python

# import libs
import csv

with open(\'$newsheet\', 'w', newline='') as outfile:
    writer = csv.writer(outfile)

    with open(\'$sheet\', 'r') as infile:
        my_reader = csv.reader(infile, delimiter=',')
        # row counter to define first line
        row_idx=0                # if first line - get index of the 3 columns needed
        datareached=0
	
        for row in my_reader:
            # read header
            if datareached == 0:
                if 'Adapter' in row:
                    writer.writerow(['AdapterRead1',row[1]])
                else:
                    writer.writerow(row)

            # if [Data] reached
            if datareached == 1:
                datareached=2

            if datareached == 2:
                if row_idx == 0:
                    sididx  = row.index('Sample_ID')
                    idxidx  = row.index('index')
                    idx2idx = row.index('index2')
                    projidx = row.index('Sample_Project')
                    row_idx = 1
                    if \'$index\' == 'dual':
                        writer.writerow(['Sample_ID','Sample_Name','Sample_Plate','Sample_Well','index','index2','Sample_Project'])
                    else:
                        writer.writerow(['Sample_ID','index','Sample_Project'])
                else:
                    currsid = row[sididx]
                    curridx = row[idxidx]
                    curridx2 = row[idx2idx]
                    currproj = row[projidx]

                    if \'$index\' == 'dual':
                        writer.writerow([currsid,currsid,'','',curridx,curridx2,currproj])
                    else:
                        writer.writerow([currsid,curridx,currproj])
            if '[Data]' in row:
                datareached=1             
	"""
}

// dragen demux
process demux {

    tag "$metaID"	
    label 'dragen'

    input:
    val newsheet from demux_sheet

    output:
    val "x" into mv_fastq

    when:
    demux = "y"
        
    """
    export LC_ALL=C        

    mkdir -p ${FQDIR}
    /opt/edico/bin/dragen --force --bcl-conversion-only=true \\
           --bcl-input-directory ${runfolder} \\
	   --output-directory ${FQDIR} \\
	   --sample-sheet ${newsheet} \\
	   --no-lane-splitting true \\
	   ${b2farg}

     """
}

process moveFastq {

    tag "${sid}_${projid}"

    input:
    val x from mv_fastq
    set sid, projid, ref from move_fastq_csv

    output:
    val "y" into run_analysis
    set sid, projid, ref into dragen_run

    when:
    demux = 'y'

    """
    mkdir -p ${OUTDIR}/${projid}
    mkdir -p ${OUTDIR}/${projid}/fastq

    mkdir -p ${OUTDIR}/${projid}/fastq/$sid

    # If there is a directory per project
    if [ -d ${FQDIR}/${projid}/ ]; then
        if [ -d ${FQDIR}/${projid}/$sid ]; then
	    mv ${FQDIR}/${projid}/$sid ${OUTDIR}/${projid}/fastq/
	# If there is a flat structure under project dir
	else
	    mv ${FQDIR}/${projid}/$sid* ${OUTDIR}/${projid}/fastq/$sid/
	fi
    # If there is a flat structure with all samples for all projects in one - create a sid folder for each sample
    else
	mv ${FQDIR}/$sid* ${OUTDIR}/${projid}/fastq/$sid/
    fi
    """

}

// Channel to start analysis if demux == 'n'
// Projects
if ( demux == 'n' ) {

   Channel
	 .fromPath(channel_sheet)
    	 .splitCsv(header:true)
    	 .map { row -> tuple( row.Sample_ID, row.Sample_Project, row.Sample_Species ) }
    	 .unique()
    	 .set{ dragen_run; }

}

// dragen run : align, vc + metrics
process dragen_align_vc {

    tag "${sid}_${projid}"
    label 'dragen' 

    input:
    val x from run_analysis
    set sid, projid, ref dragen_run

    output:
    val x into done_analyze
    val projid into multiqc_ch
    val projid  into dragen_summary    

    """
    export LC_ALL=C

    R1=\$(echo ${OUTDIR}/${projid}/fastq/${sid}/${sid}*_R1_*fastq.gz)
    R2=\$(echo ${OUTDIR}/${projid}/fastq/${sid}/${sid}*_R2_*fastq.gz)

    outdir=${OUTDIR}/${projid}/dragen/${sid}
    mkdir -p \$outdir

    /opt/edico/bin/dragen -f -r /staging/human/reference/$ref \\
        -1 \${R1} \\
    	-2 \${R2} \\
	--RGID ${projid}_${sid} \\
    	--RGSM $sid \\
	--intermediate-results-dir /staging/tmp/ \\
        --enable-map-align true \\
        --enable-map-align-output true \\
        --output-format bam \\
        --output-directory \$outdir \\
        --enable-variant-caller true \\
        --enable-sv true \\
        --output-file-prefix $sid
    
    """
}

process dragen_stats {

        tag "${projid}"

	input: 
	val projid  from dragen_summary.unique()
	
	output:
	val projid into multiqc_dragen
	val "x" into dragstats_completed
	val "x" into completed
	
	"""

	mkdir -p ${OUTDIR}/$projid/qc/dragen

	${basedir}/bin/ctg-dragen-stats-wgs -p $projid -i ${OUTDIR}/$projid/dragen/ -o ${OUTDIR}/$projid/qc/dragen/ 
	
	"""
}

process multiqc_dragen {

    tag "${projid}"

    input:
    val projid from multiqc_ch.unique()
    val x from dragstats_completed.collect()

    """
    
    cd $OUTDIR
    multiqc -f ${OUTDIR}/$projid/ --outdir ${OUTDIR}/$projid/qc/multiqc/ -n ${projid}_wgs_dragen_report.html

    mkdir -p ${CTGQC}
    mkdir -p ${CTGQC}/$projid

    cp -r ${OUTDIR}/$projid/qc ${CTGQC}/$projid/

    """

}


process multiqc_run {

	tag "$metaID"

	input:
	val x from completed.collect()

	"""

	cd $OUTDIR
	mkdir -p ${OUTDIR}/qc

	multiqc -f ${OUTDIR} $params.runfolder/ctg-interop --outdir ${OUTDIR}/qc/ -n ${metaID}_wgs_dragen_run_multiqc_report.html

	mkdir -p ${CTGQC}

	cp -r ${OUTDIR}/qc ${CTGQC}
	"""
}