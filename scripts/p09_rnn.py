

import pandas 
import os 
import pandas as pd
import numpy as np
import torch
from torch.nn import RNN, GRU, LSTM
import wandb


entity = 'your-entity-name'
project = 'your-project-name'
repo_path = 'your-repo-path'

# change working directory to the project's root
os.chdir('..')
print(os.getcwd())


# model hyperparams


# fixed params
# DATA_PATH = 'sims_out/'
PARAM_SETTINGS = range(1, 9)  # 8, starting from 1
SEEDS = range(1, 1001)   # 50, starting from 1
FIRST_TIMESTEP = 0
ALL_TIMESTEPS = [0, 1, 2, 4, 8, 12, 16]
ALL_MEASUREMENTS = ['Y1', 'Y2', 'Y3', 'Y4', 'Y5', 'Y6', 'Y7', 'Y8', 'Y9', 'Y10', 'Y11', 'Y12']
L = len(ALL_TIMESTEPS)  # number of timesteps
INPUT_SIZE = 13  # 12 measurements + 1 treatment
OUTPUT_SIZE = 12  # 12 measurements

# model hyperparams 
WHICH_RNN = 'lstm'  # 'naive' or 'gru' or 'lstm'
HIDDEN_SIZE = 50

# training hyperparams
BATCH_SIZE = 64
LR = 0.0003
N_EPOCHS = 2000

# config dictionary
# should have been done before
config = {
    'PARAM_SETTINGS': PARAM_SETTINGS,
    'SEEDS': SEEDS,
    'FIRST_TIMESTEP': FIRST_TIMESTEP,
    'ALL_TIMESTEPS': ALL_TIMESTEPS,
    'ALL_MEASUREMENTS': ALL_MEASUREMENTS,
    'L': L,
    'INPUT_SIZE': INPUT_SIZE,
    'OUTPUT_SIZE': OUTPUT_SIZE,
    # 'WHICH_RNN': WHICH_RNN,
    'HIDDEN_SIZE': HIDDEN_SIZE,
    'BATCH_SIZE': BATCH_SIZE,
    'LR': LR,
    'N_EPOCHS': N_EPOCHS,
}

# other
DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'



def load_data_setting(param_setting, seed): 
    if DEVICE == 'cpu': 
        DATA_PATH = os.path.join(repo_path, 'latent-factors/data/sims_out', "CSVs_sim_" + str(param_setting), "Y_" + str(seed) + ".csv")
    else: 
        DATA_PATH = os.path.join(repo_path, 'data/sims_out', "CSVs_sim_" + str(param_setting), "Y_" + str(seed) + ".csv")
    df = pandas.read_csv(DATA_PATH)

    return df


def dataframe_to_tensors(df):
    unique_ids = df['ID'].unique()

    # convert binary variable in TRT to integer
    df['TRT'] = df['TRT'].replace({'Placebo': 0, 'Active': 1})
    df['TRT'] = pd.to_numeric(df['TRT'])
    
    input_list, target_list, trt_list = [], [], []
    for id in unique_ids: 
        df_patient = df[df['ID'] == id]

        # inputs
        # extract TRT and first timestep as inputs
        inputs = df_patient[df_patient['Time'] == FIRST_TIMESTEP][ALL_MEASUREMENTS + ['TRT']].values
        # repeat inputs along 0-th dimension 
        # Note: this means we repeat the RNN function L-1 times, even though we use the same inputs
        inputs = inputs.repeat(L-1, axis=0)
        input_list.append(inputs)

        # targets
        # extract 2nd to last timestep as targets
        targets = df_patient[df_patient['Time'].isin(ALL_TIMESTEPS[1:])][ALL_MEASUREMENTS].values
        target_list.append(targets)

        # treatment
        trt = df_patient[df_patient['Time'] == FIRST_TIMESTEP]['TRT'].values
        trt_list.append(trt)
    
    # concatenate all inputs into a tensor, stacked across new 0 dimension; same for targets
    inputs = np.stack(input_list, axis=0)  # (200, 6, 13)
    targets = np.stack(target_list, axis=0)  # (200, 6, 12)
    trt = np.stack(trt_list, axis=0)  # (200, 1)

    # convert to torch float
    inputs = torch.tensor(inputs, dtype=torch.float32)
    targets = torch.tensor(targets, dtype=torch.float32)
    trts = torch.tensor(trt, dtype=torch.float32)

    return inputs, targets, trts

    print("here")


def load_true_ate(param_setting):
    # load true ATE
    if DEVICE == 'cpu':
        true_ate_df = pd.read_csv(os.path.join(repo_path, 'latent-factors/data/sims_out/true_ATEs.csv'))
    else: 
        true_ate_df = pd.read_csv(os.path.join(repo_path, 'data/sims_out/true_ATEs.csv'))
    true_ate = true_ate_df.iloc[param_setting-1]
    true_ate = torch.tensor(true_ate, dtype=torch.float32)

    return true_ate
     


class RNN_model(torch.nn.Module): 
    def __init__(self, which_rnn, output_size) -> None:
        super().__init__()

        # instantiate an RNN
        if which_rnn == 'naive':
            self.rnn = RNN(input_size=INPUT_SIZE, hidden_size=HIDDEN_SIZE, batch_first=True)
        elif which_rnn == 'gru':
            self.rnn = GRU(input_size=INPUT_SIZE, hidden_size=HIDDEN_SIZE, batch_first=True)
        elif which_rnn == 'lstm':
            self.rnn = LSTM(input_size=INPUT_SIZE, hidden_size=HIDDEN_SIZE, batch_first=True)
        else:
            raise ValueError("Invalid RNN type")
        
        self.fc = torch.nn.Linear(HIDDEN_SIZE * (L-1), OUTPUT_SIZE * (L-1))

    def forward(self, input):
        output, hidden = self.rnn(input)  # (batch_size, L-1, 12)
        output = torch.flatten(output, start_dim=1)
        output = self.fc(output)  
        output = output.view(-1, L-1, OUTPUT_SIZE)  # (batch_size, L-1, 12)

        return output



def run(): 
    # login wandb
    wandb.login()
    run = wandb.init(project=project, entity=entity)  # TODO anonymise

    # default config (for when not running sweep)
    sweep_config = {
        'param_setting': 1,
        'seed': 1,
        'which_rnn': 'lstm',
    }
    # set values from sweep config
    for (key, value) in wandb.config.items():
        sweep_config[key] = value
        # setattr(H, key, value)  # args is namespace object

    # add to config
    config['param_setting'] = sweep_config['param_setting']
    config['seed'] = sweep_config['seed']
    config['which_rnn'] = sweep_config['which_rnn']
    # log the config in wandb config
    wandb.config.update(config)

    print("param_setting, seed: ", sweep_config['param_setting'], sweep_config['seed'])

    df = load_data_setting(param_setting=sweep_config['param_setting'], seed=sweep_config['seed'])
    gt_ate = load_true_ate(param_setting=sweep_config['param_setting'])
    # columns: ['ID', 'Time', 'TRT', 'Y1', 'Y2', 'Y3', 'Y4', 'Y5', 'Y6', 'Y7', 'Y8', 'Y9', 'Y10', 'Y11', 'Y12']
    inputs, targets, trts = dataframe_to_tensors(df)
    if DEVICE == 'cuda':
        inputs, targets, trts = inputs.cuda(), targets.cuda(), trts.cuda()

    # normalise data
    inputs = (inputs - inputs.mean()) / inputs.std()
    targets = (targets - targets.mean()) / targets.std()

    # initialise a data loader for inputs and targets
    dataset = torch.utils.data.TensorDataset(inputs, targets, trts)

    model = RNN_model(which_rnn=sweep_config['which_rnn'], output_size=OUTPUT_SIZE).to(DEVICE)
    
    # initialise an optimiser
    optimiser = torch.optim.Adam(model.parameters(), lr=LR)

    # N_STEPS_PER_EPOCH = len(dataset)

    # a training loop
    i = 0
    for _ in range(N_EPOCHS): 
        dataloader = torch.utils.data.DataLoader(dataset, batch_size=BATCH_SIZE, shuffle=True)

        placebo_y_1_list, active_y_1_list = [], []
        for input, target, trt in dataloader:
            # forward pass
            output = model(input)
            # compute loss
            mse_loss = torch.nn.MSELoss()
            loss = mse_loss(output, target) 
            # do gradient step
            loss.backward()
            optimiser.step()
            optimiser.zero_grad()

            # wandb log; less often
            if i % 300 == 0: 
                wandb.log({'loss': loss})
            i += 1

            # Note: computed "in-distribution"
            y_1 = output[:, -1, 0]
            placebo_y_1_list.append(y_1[torch.flatten(trt == 0)])
            active_y_1_list.append(y_1[torch.flatten(trt == 1)])
        
        placebo_y_1, active_y_1 = torch.cat(placebo_y_1_list), torch.cat(active_y_1_list)
        estimated_treatment_effect = torch.mean(active_y_1) - torch.mean(placebo_y_1)

        # eval_mse = torch.square(estimated_treatment_effect - gt_ate)
        # eval_bias = estimated_treatment_effect - gt_ate
        # eval_sd = 


        # log epoch
        wandb.log({'estimated_treatment_effect': estimated_treatment_effect, 'gt_ate': gt_ate})


    print("DONE!")


if __name__ == "__main__":
    run()


    # TODO fix refactor
    # TODO all runs in parallel