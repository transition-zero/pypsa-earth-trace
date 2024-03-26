import os
import pypsa
import pandas as pd
import seaborn as sns

from helpers import get_country_networks_from_bucket


class CountryBenchmarking:

    def __init__(self, country_iso, fetch=False):

        # ---
        # Load data

        self.country_iso = country_iso
        self.historic_data = pd.read_csv('../data/ember_electricity_data.csv').query('entity_code == @country_iso')

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
        return (
            pd
            .concat([
                    self.n_unconstr.generators_t.p.filter(regex=tech).sum(axis=1).resample(resample).sum(),
                    self.n_annual_matching.generators_t.p.filter(regex=tech).sum(axis=1).resample(resample).sum(),
                ], 
                axis=1, 
                keys=['Unconstrained','Annual Matching']
            )
            .assign(historic=self.historic_data.query(" series == @tech and date.str.contains(@year) ").generation_twh.mul(1e6).values)
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
        # initiate plot
        f, ax = plt.subplots()
        # plot
        country_data.plot_total_tech_generation(tech = 'coal', year = 2019, resample='M', ax=ax)
        # formatting
        ax.set_xlabel('Month')
        ax.set_ylabel('Generation (MWh)')
        ax.set_title('Coal generation in {}'.format(iso_code))
        # save
        f.savefig('../_TRACE_outputs/generation-plots/coal_generation_{}.png'.format(iso_code), bbox_inches='tight')
        