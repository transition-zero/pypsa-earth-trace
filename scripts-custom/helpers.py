import numpy as np
import pandas as pd
import country_converter as coco

from google.cloud import storage


def get_successful_runs_from_bucket(
        bucket_name='feo-pypsa-staging'
    ) -> list:
    '''Get a list of successful runs from the bucket
    '''
    # initiate storage client
    storage_client = storage.Client()
    # get bucket
    bucket = storage_client.get_bucket(bucket_name)
    # return list
    return \
        list( 
            set([ \
                i.name.replace('results/','')[0:2] \
                    for i in bucket.list_blobs() \
                        if 'results' in i.name
                ]
        )
    )


def get_historical_data(
        path_to_data : str,
    ) -> pd.DataFrame:
    '''Get Ember's historical generation data
    '''
    historical_generation = ( 
        pd.read_csv(
            path_to_data,
        )
        .query(' Category == "Electricity generation" ')
        .query(' Subcategory == "Fuel" ')
        .query(' Unit == "TWh" ')
        .dropna(axis=0, subset=['Country code'])
        [['Area', 'Country code','Year', 'Variable', 'Unit', 'Value']]
    )

    # change iso codes from three letter to two letter
    three_letter_iso = historical_generation['Country code'].unique()
    two_letter_iso = coco.convert(three_letter_iso, to="ISO2")
    iso_mapping = {three_letter_iso[i]: two_letter_iso[i] for i in range(len(three_letter_iso))}
    historical_generation['Country code'] = historical_generation['Country code'].map(iso_mapping)

    # remap variable names
    tech_mapping = {
        'Bioenergy' : 'biomass', 
        'Coal' : 'coal', 
        'Gas' : 'gas', 
        'Hydro' : 'hydro', 
        'Nuclear' : 'nuclear', 
        'Other Fossil' : np.nan,
        'Other Renewables' : np.nan, 
        'Solar' : 'solar', 
        'Wind' : 'wind',
    }

    # overwrite existing variable names
    historical_generation['Variable'] = historical_generation['Variable'].map(tech_mapping)

    return (
        historical_generation
        .dropna(axis=0, subset=['Variable'])
        .set_index(['Country code', 'Year'])
    )