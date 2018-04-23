import pandas as pd
import numpy as np
import hashlib
from sklearn.model_selection import train_test_split
import ursgal
import pymzml
from statsmodels.compat.python import range
from statsmodels.compat.collections import OrderedDict
import os
import sys
import shutil

def main():

        file_1 = ""
        file_2 = ""
        def extract_run(files):
            for file in files:
                if file== file1:
                    run1= pymzml.run.Reader(file1, MS1_Precision =5e-6, MSn_Precision = 20e-6)
                    return run1
                if file == file2:
                    run2= pymzml.run.Reader(file2, MS1_Precision =5e-6, MSn_Precision = 20e-6)
                    return  run2

    
        def spectra_frame(runs):
            for run in runs:
                if run == run1:
                    for spectrum in run:
                        if spectrum['ms level'] ==1:
                            y = spectrum.i
                            x = spectrum.mz
                            spectra_1 = pd.DataFrame({'Run_1:M/Z':x, "Run_1:Signal_Intensity":y})
                if run == run2:          
                    for spectrum in run:
                        if spectrum['ms level'] ==1: 
                            y = spectrum.i
                            x = spectrum.mz 
                            spectra_2= pd.DataFrame({'Run_2:M/Z':x, "Run_2:Signal_Intensity":y})
                            df = pd.concat([spectra_1, spectra_2], axis = 1) 
                            return df


            def remove_nan(x):
                        return x.dropna(how = 'any')
            # add identity column to LC-MS spectra
            spectrum = []
            
            def add_identity():
                for df in  spectrum:
                    if df == spectra_1  
                                ID_1 =pd.Series([1] * len(spectra_1.index))
                                return ID_1.values
                    else:
                        if df == spectra_2  
                                ID_2 =pd.Series([1] * len(spectra_2.index))
                                return ID_2.values

            # add identity column to LC-MS/MS spectra

            result_dataframes = []
            def add_identity_2():
                for df in  result_dataframes:
                    if df == df1 
                                ID =pd.Series([1] * len(df1.index))
                                return ID.values
                    else:
                        if df == df2  
                                ID =pd.Series([1] * len(df2.index))
                                return ID.values

             # concatenate the two new dataframes together
       


            # alternative split_train_function                       
    
            def split_train_test_1(dataframe, test_ratio):
                shuffled_indices = np.random.permutation(len(dataframe))
                test_size = int(len(dataframe)*test_ratio)
                test_indices = shuffled_indices[:test_size]
                train_indices = shuffled_indices[test_size]
                return dataframe.iloc[train_indices], dataframe.iloc[test_indices]

            def get_train_test(df, y_col, x_cols, ratio):

            mask = np.random.rand(len(df)) > ratio
            df_train = df[mask]
            df_test = df[~mask]
            
            Y_train = df_train[y_col].values
            Y_test = df_test[y_col].values
            X_train = df_train[x_cols].values
            X_test = df_test[x_cols].values
            return df_train, df_test, X_train, Y_train, X_test, Y_test
        
            y_col = 
            x_cols = list(df.columns.values)
            x_cols.remove(y_col)
                
            train_test_ratio = 0.7
            df_train, df_test, X_train, Y_train, X_test, Y_test = get_train_test(df, y_col, x_cols, train_test_ratio)

            def test_set_check(identifier, test_ratio, hash):
            return hash(np.int64(identifier)).digest()[-1]< 100 *test_ratio

            def split_train_test_by_id(df, test_ratio, id_column,hash = hashlib.md5):
                ids = df.index.values
                in_test_set = ids.apply(lambda id_:test_set_check(id_, test_ratio, hash))
                return df.loc[~ in_test_set], df.loc[in_test_set]

    
    #preparation of dataframe

        df_with_id = df.reset_index()
        train_set,test_set = split_train_test_by_id(df_with_id, 0.2, "index")

    # one hot label encoder
        from sklearn.preprocessing import LabelEncoder
        encoder = LabelEncoder()
        df4_category = df5[df5['Sequence']]
        df4_category_encoded = encoder.fit_transform(df4_category)



    # Extraction of LC-MS/MS Data
                spec_files     = []

                uc = ursgal.UController(
                    params = {'database': 'Etanercept.fasta',
                            'frag_mass_tolerance': 20,
                            'frag_mass_tolerance_unit':'ppm',
                            'modifications' : [
                                'M,opt,any,Oxidation',        
                                'C,fix,any,Carbamidomethyl',  
                                '*,opt,Prot-N-term,Acetyl',],})
                results = []
                for spec_file in spec_files:
                        result = uc.search(
                            input_file = spec_file,
                            engine     = 'omssa_2_1_9',
                            force = True,)
                        results.append( result )

                        validation_engines = ['percolator_2_08', 'qvality']
                        _validated_result = uc.validate(input_file = result, engine = validation_engines,)

                        unified_file_list = []

                        unified_search_result_file = uc.search(
                            input_file = spec_files,
                            engine     = 'omssa_2_1_9',
                            force      = False
                        )
                        unified_file_list.append(unified_search_result_file)


                        merged_dataframe=pd.DataFrame()
                        extracted_dataframe = []
                        search_engine_result = pd.read_csv(result)
                        for result in search_engine_result:
                            extracted_data = search_engine_result.iloc[:, ['OMSSA:evalue', 'Exp m/z', 'OMSSA:pvalue', 'Retention Time(s)']]
                            extracted_dataframe.append(extracted_data)
                            merged_dataframe = pd.concat(extracted_dataframe, axis = 1)
                            _corrected_merged_frame = merged_dataframe.fillna(value = 0)


                            uc.visualize(input_files    = unified_file_list, engine= 'venndiagram',)
                            return

if __name__ == '__main__':
    main()   
        
 