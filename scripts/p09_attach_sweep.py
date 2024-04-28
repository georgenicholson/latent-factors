import wandb
import argparse

from p09_rnn import run

parser = argparse.ArgumentParser(description='Attaching to a sweep.')
parser.add_argument('--sweep_id', type=str, default="4nvmgvhv", metavar='N', help="ID of sweep to attach to.")
H, unknown = parser.parse_known_args()

wandb.agent("fabianfalck" + '/' + "jbi_second" + '/' + H.sweep_id, function=run)