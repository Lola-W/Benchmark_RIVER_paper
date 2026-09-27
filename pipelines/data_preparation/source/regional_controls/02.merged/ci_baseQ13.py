import os
import sys
import subprocess



python_script = "./work/regional_controls/02.merged/compute_maf_binomial_mnv.py"
header = "\t".join(["#CHROM", "POS", "REF", "ALT", \
                    "REF_COUNT", "ALT_COUNT", \
                    "MAF", "LOWER_CI", "UPPER_CI", "IS_GREATER"])


def main(argv):
    if len(argv) != 6:
        sys.stderr.write("usage: " + argv[0] + " <fasta> <input.bam> <input.vcf> <samtools_output_dir> <output.vcf>\n")
        sys.exit(2)
    fasta = argv[1]
    input_bam = argv[2]
    input_vcf = argv[3]
    samtools_output_dir = argv[4]
    if not samtools_output_dir.endswith("/"):
        samtools_output_dir += "/"
    if not os.path.exists(samtools_output_dir):
        os.makedirs(samtools_output_dir)
    output_vcf = argv[5]
    
    wf = open(output_vcf, "w")
    rf = open(input_vcf, "r")
    wf.write(header + "\n")
  
    for line in rf:
        if line.startswith("#"):
            continue
        line = line.rstrip().split("\t")
        print(line)
        chrom = line[0]
        pos = line[1]
        ref = line[3]
        alt = line[4]
        start_pos = str(int(pos) - 5)
        # if len(pos) > len(alt):
        #     end_pos = str(int(pos) + len(ref) + 4)
        # else:
        #     end_pos = str(int(pos) + 5)
        start_pos = int(pos) - 5
        end_pos = int(pos) + max(len(ref), len(alt)) - 1 + 5
        positions = f"{chrom}:{start_pos}-{end_pos}"
        # positions = chrom + ":" + start_pos + "-" + end_pos

        samtools_outfile = samtools_output_dir +  positions + ".txt"
        samtools_command = "samtools mpileup -r " + positions + \
                           " -f " + fasta + " -Q13 -q0 -AB --ff SECONDARY,QCFAIL,DUP,SUPPLEMENTARY -d1000000 " \
                           + input_bam + " > " + samtools_outfile
        subprocess.call(samtools_command, shell=True)
        if os.path.exists(samtools_outfile):
            python_command = ["python", python_script, samtools_outfile, chrom + ":" + pos, ref, alt]
            output,error = subprocess.Popen(python_command,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
            output = output.decode()
            print(output)
            try:
                num_ref, num_alt, maf, lower_ci, upper_ci, is_greater = output.rstrip().split(" ")
                wf.write("\t".join([chrom, pos, ref, alt, num_ref, num_alt, maf, lower_ci, upper_ci, is_greater]) + "\n")
            except ValueError:
                # If output format is unexpected, write what we have
                wf.write("\t".join([chrom, pos, ref, alt, output.rstrip()]) + "\n")
        else:
            wf.write("\t".join([chrom, pos, ref, alt, "nan", "nan", "nan", "nan", "nan", "nan"]) + "\n")

    wf.close()
    rf.close()


main(sys.argv)      

