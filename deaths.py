import pandas as pd
import json
import numpy as np

# df = pd.read_json('deaths.json')

with open('deaths.json') as file:
    file_data = file.read()
data = json.loads(file_data)

index = data['labels']
index

len(index)

columns = {'index': index}

for column in data['datasets']:
    columns[column['label']] = column['data']

columns.keys()

df = pd.DataFrame(columns, dtype=int)
df.plot()
log_df = df.rolling(7).mean().apply(np.log10)

log_df.sub(log_df['80+ Jahre'], axis=0).plot()


# Need nowcasting of likely cases that will still be gemeldet
# Need death according to date of death Sterbedatum vs. Meldedatum vs. in RKI stats

# Where can I find deaths according to date of death?
# Then check historically how many cases were recorded in the days thereafter to nowcast.
