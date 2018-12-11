import csv
import argparse

parser = argparse.ArgumentParser(description='Extract MLST Profile')
parser.add_argument('--profile', '-p', help='mlst profile', required=True)
parser.add_argument('--blast', '-b', help='mlst profile', required=True)
args = parser.parse_args()


def load_profiles(profile_file):
    """Load mlst profile file, format profiles to dictionaries."""
    def identify_columns(header):
        return [header[i] for i in range(1, len(header) - 1)]

    def format_profiles(line, header):
        formated_profile = {}
        line = line[1:-1]
        for i in range(len(line)):
            formated_profile[header[i]] = line[i]
        return formated_profile

    with open(profile_file) as csvfile:
        rdr = csv.reader(csvfile, delimiter='\t')
        lines = [line for line in rdr]
        header = identify_columns(lines[0])
        profiles = {line[0]: format_profiles(line, header) for line in lines[1:]}

    return profiles


def load_blast_result(blast_result):
    """Load blast result file, extract gene alleles."""
    with open(blast_result) as csvfile:
        rdr = csv.reader(csvfile, delimiter='\t')
        lines = [line for line in rdr]
    mlst_profile = {line[0].split('_')[0]: line[0].split('_')[1] for line in lines}

    return mlst_profile


def output(profile_file, blast_result):
    """Find strain type from mlst profile."""
    all_profiles = load_profiles(profile_file)
    mlst_profile = load_blast_result(blast_result)

    for st_type, profile in all_profiles.items():
        if mlst_profile == profile:
            print('Strain Type:', st_type, sep=' ')


if __name__ == '__main__':
    output(args.profile, args.blast)
