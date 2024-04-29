import wandb
import traceback


from p09_rnn import run

entity = 'your-entity-name'
project = 'your-project-name'


sweep_config = {
    'name': 'my-sweep',
    'method': 'grid',  # or 'grid' or 'bayes'

    'parameters': {
        'param_setting': {
            'values': list(range(1, 9)),  # ,
        },
        'seed': {
            'values': list(range(1, 1001)),
        },
        'which_rnn': {
            'values': ['naive', 'gru', 'lstm'],
        },
    }
}



sweep_id = wandb.sweep(sweep_config, project=project, entity=entity)
wandb.agent(sweep_id, function=run)

