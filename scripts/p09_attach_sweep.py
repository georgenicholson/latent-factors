import wandb
import argparse

from p09_rnn import run

entity = 'your-entity-name'
project_name = 'your-project-name'

parser = argparse.ArgumentParser(description='Attaching to a sweep.')
parser.add_argument('--sweep_id', type=str, default="4nvmgvhv", metavar='N', help="ID of sweep to attach to.")
H, unknown = parser.parse_known_args()

wandb.agent(entity + '/' + project_name + '/' + H.sweep_id, function=run)