import numpy as np
import pandas as pd 
import subprocess
import yaml
import secrets

user_df = pd.read_csv('credentials.csv', header=None, usecols=[0,1,2], names=['username','password','api_key'])
user_df['username'] = user_df['username'].map(lambda x: x.lower())

# user_df['api_key'] = user_df['username'].map(lambda x: secrets.token_urlsafe(16))
# user_df['api_key'] = 'c7235440-f6b3-4f31-b0b0-ff323f99a9a7'
user_df['workspaces'] = '' # user_df['username'].map(lambda x: ['cancer_stage_1', 'cancer_stage_2', 'cancer_stage_3'] if x not in ['davidkartchner','alexisnunn','cassiemitchell','haydnturner','dongyuzhang'] else None)

# Get users that have been recently added
old_users = pd.DataFrame(yaml.safe_load(open('.users.yaml','r')))['username'].tolist()
new_user_df = user_df[user_df.username.map(lambda x: x not in old_users)]
new_user_df['hashed_password'] = new_user_df.password.map(lambda x: subprocess.check_output(f'htpasswd -nbB "" {x}'.split()).decode("utf-8").strip().strip('"').strip(':'))

records = new_user_df[['username','hashed_password', 'api_key','workspaces']].to_dict(orient='records')



yaml.dump(records, open('.new_users.yaml', 'w'))
# print(user_df[['username','hashed_password']].to_dict(orient='records'))