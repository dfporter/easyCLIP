from typing import List, Mapping, Union
import pandas


def get_name_translator(
    fname: str=f'/Users/dp/pma/RBP missense mutations/Sequencing_schemes_and_results.xlsx',
    sheet_name='Proteins'):

    to_name = {'100P': 'PCBP1 L100P', '100Q': 'PCBP1 L100Q',
              'PCBP1 100P': 'PCBP1 L100P', 'PCBP1 100Q': 'PCBP1 L100Q',
              'RQCD1': 'CNOT9', 'RQCD1 MUT': 'CNOT9 P131L', 'NUFIP': 'NUFIP1',
              'NUFIP1 MUT': 'NUFIP1 R475W',
              'NUFIP MUT': 'NUFIP1 R475W', 'YTHDC2 635K': 'YTHDC2 E635K',
              'YTHDC2 185K': 'YTHDC2 E185K'}

    for name in ['SF3B1 E902K', 'SF3B1 R625C', 'SF3B1 R625H', 'SF3B1 K700E',
        'YTHDC2 E185K', 'YTHDC2 E635K']:
        to_name[name] = name

    def determine_basename(name):
        if type(name) != type(''):
            return ''
        if name == 'hnRNP C' or name == 'hnRNP D':
            return name
        else:
            return name.split(' ')[0]
    
    def pair(name):
        if type(name) != type(''):
            return 
        
        if ' ' in name and (determine_basename(name) != name):
            to_name[determine_basename(name) + ' MUT'] = name
            to_name[determine_basename(name) + ' mut'] = name
            to_name[determine_basename(name) + ' Mut'] = name
        else:
            to_name[name] = name

    df = pandas.read_excel(fname, sheet_name=sheet_name)
    [pair(x) for x in df['Proteins'].tolist()]
    [pair(x) for x in df['Mutations'].tolist()]

    return to_name

to_name = get_name_translator()

def translate_to_proper_names(name_list: List[str]) -> List[str]:
    return [to_name.get(name, '') for name in name_list]

def get_mutant_for_a_name(name: str) -> str:
    if ' ' in name:
        return to_name[name]
    return to_name[name + ' MUT']    

def get_mutant(name_or_name_list: Union[List[str], str]) -> Union[List[str], str]:
    print('get_mut(): ', name_or_name_list)
    if type(name_or_name_list) == type(''):
        return get_mutant_for_a_name(name_or_name_list)

    return [get_mutant_for_a_name(name) for name in name_or_name_list]
        