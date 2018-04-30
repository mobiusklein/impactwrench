import ursgal
import os
import os
import click
import os
import glob
from os import listdir


input_file =click.prompt('Please enter a valid file path to a fasta_file', type=str)

params = {
    'enzyme': 'trypsin',
    'decoy_generation_mode': 'reverse_protein',
}

uc = ursgal.UController(
    params=params)

new_target_decoy_db_name = uc.execute_misc_engine(
    input_file=input_file,
    engine='generate_target_decoy_1_0_0',
    output_file_name='Etanercept_target_decoy.fasta',)