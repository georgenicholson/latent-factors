import wandb
import pickle
import os
import torch

# Your sweep ID !!!!!
sweep_id = 'zpkqra53'
SAVE = False

project = 'your-project-name'
entity = 'your-entity-name' # Typically your username or team name in wandb


PARAM_SETTINGS = range(1, 9)  # 8, starting from 1

if SAVE: 

    # create sweep folder in ..data/sweep_id
    os.makedirs(f'results/{sweep_id}', exist_ok=True)

    # Initialize wandb API
    api = wandb.Api()

    # Fetch all runs in the sweep
    sweep_runs = api.sweep(f"{entity}/{project}/{sweep_id}").runs

    # Dictionary to hold metrics and configurations
    # runs_data = {} 

    # Iterate through all runs and retrieve metrics
    for i, run in enumerate(sweep_runs):
        if i % 100 == 0: 
            print("run i: ", i)
        # Fetch configuration for the run
        config = run.config
        run_time = run.summary['_runtime']

        # Fetch all logged metrics for a run
        # if config['param_setting'] != 8:   # some issue with 8
        metrics = {metric: run.history(keys=[metric]) for metric in ['estimated_treatment_effect', 'gt_ate']}

        # slow V2
        # # Fetch the full history (this contains all metrics logged during the run)
        # history_df = run.history()
        
        # # Extract specific metrics from history_df as needed. For example:
        # # If you are interested in a metric named 'accuracy':
        # # metrics = history_df['accuracy'].tolist() if 'accuracy' in history_df.columns else None
        
        # # For demonstration, let's assume we collect all metrics directly
        # metrics = history_df.to_dict('list')  # Converts the dataframe to a dictionary of lists


        # Store metrics and configuration
        # runs_data[run.id] = {'metrics': metrics, 'config': config}

        # pickle the metrics and config 
        run_id = run.id
        with open(f'results/' + sweep_id + '/' + run_id + '.pkl', 'wb') as f:
            pickle.dump({'metrics': metrics, 'config': config, 'run_time': run_time}, f)



# load the metrics and config from results/ directory
naive_list, gru_list, lstm_list = [], [], []
naive_dict, gru_dict, lstm_dict = {param_setting: [] for param_setting in PARAM_SETTINGS}, {param_setting: [] for param_setting in PARAM_SETTINGS}, {param_setting: [] for param_setting in PARAM_SETTINGS}
naive_time_list, gru_time_list, lstm_time_list = [], [], []
for file in os.listdir('results/' + sweep_id):
    if file.endswith('.pkl'):
        with open('results/' + sweep_id + '/' + file, 'rb') as f:
            res_dict = pickle.load(f)
            estimated_treatment_effect = res_dict['metrics']['estimated_treatment_effect']['estimated_treatment_effect'].values[-1]
            try: 
                gt_ate = res_dict['metrics']['gt_ate']['gt_ate'].values[-1]
            except:
                 gt_ate = 123456789

            if res_dict['config']['which_rnn'] == 'naive':
                naive_dict[res_dict['config']['param_setting']].append([estimated_treatment_effect, gt_ate])
                naive_time_list.append(res_dict['run_time'])
            elif res_dict['config']['which_rnn'] == 'gru':
                gru_dict[res_dict['config']['param_setting']].append([estimated_treatment_effect, gt_ate])
                gru_time_list.append(res_dict['run_time'])
            elif res_dict['config']['which_rnn'] == 'lstm':
                lstm_dict[res_dict['config']['param_setting']].append([estimated_treatment_effect, gt_ate])
                lstm_time_list.append(res_dict['run_time'])

            # all together
            res = [estimated_treatment_effect, gt_ate]
            if res_dict['config']['which_rnn'] == 'naive':
                naive_list.append(res)
            elif res_dict['config']['which_rnn'] == 'gru':
                gru_list.append(res)
            elif res_dict['config']['which_rnn'] == 'lstm':
                lstm_list.append(res)

            # print(res_dict)
            # print("metrics: ", data['metrics'])
            # print("config: ", data['config'])

# convert lists to tensors
naive_tensor = torch.tensor(naive_list)
gru_tensor = torch.tensor(gru_list)
lstm_tensor = torch.tensor(lstm_list)

names = ['naive', 'gru', 'lstm']

for tensor, name in zip([naive_tensor, gru_tensor, lstm_tensor], names):
    mse = torch.mean((tensor[:, 0] - tensor[:, 1])**2)
    bias = torch.mean(tensor[:, 0] - tensor[:, 1])
    std = torch.std(tensor[:, 0] - tensor[:, 1])

    print(f"RNN type: {name}, MSE: {mse}, Bias: {bias}, Std: {std}")


# save to csv files
for param_setting in naive_dict.keys(): 
    with open(f'results/' + sweep_id + '_' + 'naive_' + str(param_setting) + '.csv', 'w') as f:
        f.write('estimated_treatment_effect,gt_ate\n')
        for res in naive_dict[param_setting]: 
                f.write(str(res[0]) + ',' + str(res[1]) + '\n')

for param_setting in gru_dict.keys():
    with open(f'results/' + sweep_id + '_' + 'gru_' + str(param_setting) + '.csv', 'w') as f:
        f.write('estimated_treatment_effect,gt_ate\n')
        for res in gru_dict[param_setting]:
                f.write(str(res[0]) + ',' + str(res[1]) + '\n')

for param_setting in lstm_dict.keys():
    with open(f'results/' + sweep_id + '_' + 'lstm_' + str(param_setting) + '.csv', 'w') as f:
        f.write('estimated_treatment_effect,gt_ate\n')
        for res in lstm_dict[param_setting]:
                f.write(str(res[0]) + ',' + str(res[1]) + '\n')

# save run times
with open(f'results/' + sweep_id + '_' + 'naive_run_times.csv', 'w') as f:
    f.write('run_time\n')
    for res in naive_time_list: 
            f.write(str(res) + '\n')
        
with open(f'results/' + sweep_id + '_' + 'gru_run_times.csv', 'w') as f:
    f.write('run_time\n')
    for res in gru_time_list: 
        f.write(str(res) + '\n')

with open(f'results/' + sweep_id + '_' + 'lstm_run_times.csv', 'w') as f:
    f.write('run_time\n')
    for res in lstm_time_list: 
        f.write(str(res) + '\n')





print("fully done.")