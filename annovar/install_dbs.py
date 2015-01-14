'''
Methods for installing ANNOVAR databases.
'''

import os
import subprocess
import argparse

# ANNOVAR databases to install.
ANNOVAR_DATABASES = [
    # -- Genes --

    # UCSC gene annotation.
    { 'name': 'refGene', 'host': 'ucsc' },

    # Gencode.
    { 'name': 'wgEncodeGencodeCompV19', 'host': 'ucsc' },

    # -- Regions --

    # Segmental duplication areas from UCSC.
    { 'name': 'genomicSuperDups', 'host': 'ucsc' },
    # Conservation.
    { 'name': 'phastConsElements46way', 'host': 'ucsc' },

    # -- Filtering databases --

    # NOTE: these next three entries are excluded by default because they are very large.

    # Most recent 1000 genomes.
    { 'name': '1000g2014oct', 'host': 'annovar' },

    # Most recent dbSNP + non-flagged as well.
    { 'name': 'snp138', 'host': 'annovar' },
    { 'name': 'snp138NonFlagged', 'host': 'annovar' },

    # COSMIC cancer db.
    { 'name': 'cosmic70', 'host': 'annovar' },

    # NCI60.
    { 'name': 'nci60', 'host': 'annovar' },

    # ClinVar.
    { 'name': 'clinvar_20140929', 'host': 'annovar' },

    # NHLBI exome sequencing project.
    { 'name': 'esp6500si_all', 'host': 'annovar' },

    # CADD, > 1% (3GB) and > 10% (33GB)
    { 'name': 'caddgt10', 'host': 'annovar' },
    #{ 'name': 'caddgt20', 'host': 'annovar' },

    # ExAC exome sequencing.
    { 'name': 'exac02', 'host': 'annovar' }
]

def install_annovar_dbs(annovar_dir='.', build='hg19', dest_dir='humandb'):
    '''
    Installs all ANNOVAR databases.
    '''

    for entry in ANNOVAR_DATABASES:
        if isinstance( entry, dict ):
            database = entry['name']
            host = entry['host']
        else:
            database = entry
            host = 'ucsc'

        install_annovar_db(database=database, annovar_dir=annovar_dir, build=build, host=host, dest_dir=dest_dir)

def install_annovar_db(database, annovar_dir='.', build='hg19', host='annovar', dest_dir='humandb'):
    '''
    Install an ANNOVAR database.
    '''
    annovar = os.path.join(annovar_dir, 'annotate_variation.pl')
    subprocess.call('%s -build %s -downdb %s -webfrom %s %s' % (annovar, build, database, host, dest_dir), shell=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Install ANNOVAR databases.')
    parser.add_argument('--database', help='Database to install')
    args = parser.parse_args()

    # Bare-bones database installation.
    if args.database:
        install_annovar_db(database=args.database)
    else:
        install_annovar_dbs()




