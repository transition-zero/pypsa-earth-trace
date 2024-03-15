import pandas as pd

def get_historical_data(
        path_to_data : str,
    ) -> pd.DataFrame:

    historical_generation = ( 
        pd.read_csv(
            path_to_data, 
            index_col = ['Country code','Year']
        )
        .query(' Category == "Electricity generation" ')
        .query(' Subcategory == "Fuel" ')
        .query(' Unit == "TWh" ')
        [['Area', 'Variable', 'Unit', 'Value']]
    )

    return historical_generation
