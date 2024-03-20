import numpy as np
import pandas as pd
import country_converter as coco

from google.cloud import storage


def get_successful_runs_from_bucket(
        bucket_name='feo-pypsa-staging'
    ) -> list:
    '''Returns a list of countries (two-letter iso codes) for which we have
    managed successfully to solve PyPSA optimisation models.
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


def download_all_prepared_networks(
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
        field_to_get : str = 'Electricity generation',
    ) -> pd.DataFrame:
    '''Get Ember's historical generation data
    '''
    historical_generation = ( 
        pd.read_csv(
            path_to_data,
        )
        .query(f' Category == "{field_to_get}" ')
        .query(' Subcategory == "Fuel" ')
        #.query(' Unit == "TWh" ')
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


def make_tracker_sheet(
        save=False,
        base_year=2019,
    ):

    historical = get_historical_data(
        path_to_data='../data/ember_electricity_data.csv',
        field_to_get='Power sector emissions'
    )

    # get iso codes for dataframe
    iso_codes = historical.index.get_level_values(0).unique()

    # make dataframe
    trace_tracker = pd.DataFrame({
        'iso' : iso_codes,
        'country' : coco.convert(names=iso_codes, to='name_short')
    })

    # append emissions
    trace_tracker[f'emissions_mtco2_{base_year}'] = \
        trace_tracker['iso'].map(
            historical
            .xs(base_year, level=1)
            .reset_index()
            .groupby(by=['Country code'])
            .sum(numeric_only=True)
            .Value
            .to_dict()
    )

    # add emissions share
    trace_tracker['emissions_share_%'] = trace_tracker[f'emissions_mtco2_{base_year}'] / trace_tracker[f'emissions_mtco2_{base_year}'].sum() * 100

    # append coal generation
    coal_generation = (
        get_historical_data(path_to_data='../data/ember_electricity_data.csv')
        .query('Unit == "TWh"')
        .query('Variable == "coal"')
        .xs(base_year, level=1)
        .Value
        .to_dict()
    )

    trace_tracker[f'coal_generation_TWh_{base_year}'] = trace_tracker['iso'].map(coal_generation).fillna(0)

    # append trace estimates
    trace_estimates = (
        pd
        .read_csv('../data/power/electricity-generation_emissions-sources.csv')
        .dropna(subset=['activity', 'activity_units'], axis=0)
        #.query(' source_type.str.contains("coal") ')
        .query(' source_type == "coal" ')
        .query(f' start_time.str.contains("{base_year}") ')
        .groupby(['iso3_country'], as_index=False)
        .sum(numeric_only=True)
    )

    trace_estimates = (
        trace_estimates
        .assign(iso = coco.convert(trace_estimates['iso3_country'], to='iso2'))
        [['iso', 'activity']]
        .set_index('iso')
        .activity
        .divide(1e6) # convert to TWh
        .to_dict()
    )

    trace_tracker[f'coal_generation_trace_current_TWh_{base_year}'] = trace_tracker['iso'].map(trace_estimates).fillna(0)

    # to csv
    if save:
        trace_tracker.to_csv('../TRACE-BENCHMARKS.csv', index=False)

    return trace_tracker