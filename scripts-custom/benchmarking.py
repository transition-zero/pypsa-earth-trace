import os
import pypsa
import pandas as pd
import matplotlib.pyplot as plt

import seaborn as sns
sns.set_theme(style="darkgrid", palette="pastel")

from helpers import (
    get_modelling_progress,
    get_country_networks_from_bucket
)


class CountryBenchmarking:

    def __init__(self, country_iso, fetch=False):

        # ---
        # Load data

        self.country_iso = country_iso

        # get historic generation data
        self.historic_generation = (
            pd.
            read_csv('../data/ember-electricity-generation-monthly.csv')
            .query('entity_code == @country_iso')
        )

        # get historic demand data
        self.historic_demand = (
            pd.
            read_csv('../data/ember-electricity-demand-monthly.csv')
            .query('entity_code == @country_iso')
        )

        # get historic capacity data
        self.historic_capacity = (
            pd.
            read_csv('../data/ember-electricity-capacity-yearly.csv')
            .query('entity_code == @country_iso')
        )

        # get latest data from bucket
        if fetch:
            get_country_networks_from_bucket(bucket_name='feo-pypsa-staging', country_iso=country_iso)

        self.n_unconstr = pypsa.Network(f'../feo-pypsa-staging/results/{country_iso}/elec_s_10_ec_lv1.00_1H_trace.nc')
        self.n_annual_matching = pypsa.Network(f'../feo-pypsa-staging/results/{country_iso}/elec_s_10_ec_lv1.00_1H-constr_trace.nc')


    def get_total_tech_generation(
            self, tech, year, resample='M',
        ):
        '''Get sum of generation for a given technology at a specified temporal resolution (e.g., 'D', 'M', 'Y')
        '''
        year = str(year)

        if tech == 'gas':
            tech_pypsa = 'OCGT|CCGT'
        else:
            tech_pypsa = tech
            
        return (
            pd
            .concat([
                    self.n_unconstr.generators_t.p.filter(regex=tech_pypsa).sum(axis=1).resample(resample).sum(),
                    self.n_annual_matching.generators_t.p.filter(regex=tech_pypsa).sum(axis=1).resample(resample).sum(),
                ], 
                axis=1, 
                keys=['Unconstrained','Annual Matching']
            )
            .assign(historic = \
                    self
                    .historic_generation
                    .query(" series == @tech and date.str.contains(@year) ")
                    .generation_twh
                    .mul(1e6)
                    .values
            )
            .reset_index()
            .melt(id_vars='snapshot', var_name='Constraint', value_name='Generation (MWh)')
        )
    
    
    def plot_total_tech_generation(self, tech, year, resample='M', ax=None):
        '''Plot total generation for a given technology
        '''
        return (
            sns.lineplot(
                data=self.get_total_tech_generation(tech, year, resample), 
                x='snapshot', 
                y='Generation (MWh)', 
                hue='Constraint',
                ax=ax,
            )
        )
    

    def get_capacity(self, year = 2019):
        '''Returns a pandas dataframe containing the historic and modelled capacity for a given year.
        '''
        return pd.concat([
            (self
                .historic_capacity
                .query("date == @year")
                .pivot_table(columns='series', values='capacity_gw', aggfunc='max')
                .set_index( pd.Index( ['Ember'] ) )
                .round(1)
            ),
            (self
                .n_unconstr
                .generators
                .groupby(by='carrier')
                .sum(numeric_only=True)
                .p_nom
                .reset_index()
                .pipe(lambda x: x.assign(
                    carrier = x.carrier.map({
                        'CCGT' : 'gas',
                        'OCGT' : 'gas',
                        'offwind-ac' : 'wind',
                        'offwind-dc' : 'wind',
                        'onwind' : 'wind',
                        'ror' : 'hydro',
                        'solar' : 'solar',
                        'nuclear' : 'nuclear',
                        'hydro' : 'hydro',
                        'coal' : 'coal',
                        'oil' : 'oil',
                        'geothermal' : 'geothermal',
                        'biomass' : 'biomass',
                    })
                    )
                )
                .dropna(axis=0)
                .pivot_table(columns='carrier', values='p_nom', aggfunc='sum')
                .mul(1e-3)
                .round(1)
                .set_index( pd.Index( ['PyPSA-Earth'] ) )
                )
        ])
    

    def plot_total_tech_capacity(self, year = 2019):
        '''Plots the total capacity for each technology in a given year.
        '''
        f, ax = plt.subplots()
        
        (self
         .get_capacity(year)
         .plot(
             kind='bar', 
             stacked=True, 
             title=f'Capacity in {self.country_iso} ({year})', 
             ylabel='Capacity (GW)',
             ax=ax,
        ))

        return f, ax

    
    def get_demand(self, year):
        year = str(year)
        return (
            pd.DataFrame({
                'Ember' : [self.historic_demand.query("date.str.contains(@year)").groupby(by='entity_code').sum(numeric_only=True).demand_twh.mul(1e6).values[0]],
                'PyPSA-Earth' : [self.n_unconstr.loads_t.p_set.sum().sum()],
            })
            .melt()
        )

    def plot_demand(self, year):
        '''Benchmarks modelled demand against historical data
        '''
        demand_data = self.get_demand(year)

        f, ax = plt.subplots()

        ax = sns.barplot(
            data=demand_data, 
            x='variable', 
            y='value',
            ax=ax
        )

        # calculate pctg change for annotation
        x = demand_data.query( " variable == 'Ember' ").value.iloc[0]
        y = demand_data.query( " variable == 'PyPSA-Earth' ").value.iloc[0]
        pctg_change = round( (y - x)/x * 100, 1)

        if pctg_change >= 0:
            symb = '+'
        else:
            symb = '-'

        # Annotate each bar with its value
        for i, p in enumerate(ax.patches):
            if i > 0:
                ax.annotate(
                    symb + str(pctg_change) + '%',
                    (p.get_x() + p.get_width() / 2., p.get_height()*0.9 ),
                    ha='center', 
                    va='bottom'
                )

        ax.set_ylabel('Electricity demand (MWh)')
        ax.set_xlabel('')
        ax.set_title(f'Electricity demand in {self.country_iso} ({year})')

        return f, ax
    

if __name__ == "__main__":

    import warnings
    warnings.filterwarnings("ignore")
    
    import matplotlib.pyplot as plt

    failed_iso = []
    for cfig in os.listdir('../country_configs'):
        # get iso code
        iso_code = cfig.split('.')[1]
        # try to fetch networks from bucket
        try:
            country_data = CountryBenchmarking(iso_code, fetch=False)
        except:
            failed_iso.append(iso_code)
            continue
        
        # ---
        # BENCHMARK GENERATION

        for tech in ['coal', 'solar', 'gas', 'nuclear']:
            try:
                # initiate plot
                f, ax = plt.subplots()
                # plot
                country_data.plot_total_tech_generation(tech = tech, year = 2019, resample='M', ax=ax)
                # formatting
                ax.set_xlabel('Month')
                ax.set_ylabel('Generation (MWh)')
                ax.set_title(f'{tech.title()} generation in {iso_code}')
                # save
                f.savefig(f'../_TRACE_outputs/generation-plots/{iso_code}_{tech}_generation.png', bbox_inches='tight')
            except:
                pass
                continue

        # ---
        # BENCHMARK DEMAND
        f, ax = country_data.plot_demand(year = 2019)
        f.savefig(f'../_TRACE_outputs/demand-plots/{iso_code}_demand_2019.png', bbox_inches='tight')

        # ---
        # BENCHMARK CAPACITY
        f, ax = country_data.plot_total_tech_capacity(year = 2019)
        f.savefig(f'../_TRACE_outputs/capacity-plots/{iso_code}_capacity_2019.png', bbox_inches='tight')
    
    # ---
    # OUTPUT MODELLING PROGRESS
    
    get_modelling_progress()