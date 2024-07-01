# %%
import random

import ase
import ase.io
import ase.io.cif
import ase.io.vasp
import ase.io.lammpsdata

from io import StringIO

import argparse

try:
    import sqsgenerator
    from sqsgenerator import from_ase_atoms
    has_sqsgenerator = True

except ImportError:
    has_sqsgenerator = False
    ## We dont need to use the bool I guess



def build_reference_system(x_dim, y_dim, z_dim):
    #### ORIGINAL CIF FILE ####
    # data_global
    # _chemical_name_mineral 'Microcline'
    # loop_
    # _publ_author_name
    # 'Allan D R'
    # 'Angel R J'
    # _journal_name_full 'European Journal of Mineralogy'
    # _journal_volume 9 
    # _journal_year 1997
    # _journal_page_first 263
    # _journal_page_last 275
    # _publ_section_title
    # ;
    #  A high-pressure structural study of microcline (KAlSi3O8) to 7 GPa
    #  P = 0.0 GPa
    # ;
    # _database_code_amcsd 0006649
    # _chemical_formula_sum '(K.986 Na.014) (Al1.03 Si2.97) O8'
    # _cell_length_a 8.5733
    # _cell_length_b 12.9375
    # _cell_length_c 7.2075
    # _cell_angle_alpha 90.530
    # _cell_angle_beta 115.972
    # _cell_angle_gamma 87.968
    # _cell_volume 718.229
    # _exptl_crystal_density_diffrn      2.572
    # _symmetry_space_group_name_H-M 'C -1'
    # loop_
    # _space_group_symop_operation_xyz
    #   'x,y,z'
    #   '1/2+x,1/2+y,z'
    #   '-x,-y,-z'
    #   '1/2-x,1/2-y,-z'
    # loop_
    # _atom_site_label
    # _atom_site_fract_x
    # _atom_site_fract_y
    # _atom_site_fract_z
    # _atom_site_occupancy
    # _atom_site_U_iso_or_equiv
    # KM   0.28410  -0.00710   0.13670   0.98600   0.02495
    # NaM   0.28410  -0.00710   0.13670   0.01400   0.02495
    # AlT1o   0.01020   0.18810   0.21800   0.94400   0.01646
    # SiT1o   0.01020   0.18810   0.21800   0.05600   0.01646
    # AlT1m   0.00880   0.81820   0.23150   0.06800   0.01646
    # SiT1m   0.00880   0.81820   0.23150   0.93200   0.01646
    # AlT2o   0.71050   0.12060   0.34050   0.00900   0.01646
    # SiT2o   0.71050   0.12060   0.34050   0.99100   0.01646
    # AlT2m   0.70690   0.88640   0.34970   0.00900   0.01697
    # SiT2m   0.70690   0.88640   0.34970   0.99100   0.01697
    # Oa1   0.00170   0.14060  -0.01570   1.00000   0.02191
    # Oa2   0.63960   0.00500   0.28590   1.00000   0.01912
    # Obo   0.82110   0.14990   0.22050   1.00000   0.02482
    # Obm   0.83240   0.85580   0.23800   1.00000   0.02634
    # Oco   0.03470   0.32140   0.25300   1.00000   0.02254
    # Ocm   0.03840   0.69120   0.26840   1.00000   0.02090
    # Odo   0.18840   0.12530   0.40450   1.00000   0.02014
    # Odm   0.17680   0.87000   0.41020   1.00000   0.02115

    ##### BUT .cif in ase doesnt like C-1 space group, lets work in lammpsdata ####

    Microcline_with_Al_Order_Logic_Embedded_File= StringIO("""
    # Microcline - (K.986 Na.014) (Al1.03 Si2.97) O8
    
            52  atoms
            6  atom types
    
        0.000000000000       8.573300000000  xlo xhi
        0.000000000000      12.929364636010  ylo yhi
        0.000000000000       6.479443178941  zlo zhi
        0.458733429218      -3.156393886927       0.045276610361  xy xz yz
    
    Masses
    
                1   15.99900000             # O
                2   39.09830000             # K
                3   276.00000000            # Mt - PLACEOHOLDER FOR T1O
                4   281.00000000            # Ds - PLACEOHOLDER FOR T1M
                5   280.00000000            # Rg - PLACEOHOLDER FOR T2O
                6   285.17000000            # Cn - PLACEOHOLDER FOR T2M
    
    Atoms # atomic
    
            1    1     -3.027765972754       1.862434435402       6.377715921031
            2    1      4.583363334874       0.077591406082       1.852472804859
            3    1      6.412315918972       1.948095251523       1.428717220956
            4    1      6.777777243636      11.075726088764       1.542107476588
            5    1     -0.353637219242       4.166952776435       1.639299124272
            6    1     -0.200884852976       8.948929078631       1.739082549228
            7    1      0.395927691419       1.638363777783       2.620934765882
            8    1      0.620104751002      11.267119698899       2.657867592002
            9    2      2.459671907528      12.843755459731       0.885739882561
            10    3     -0.514358449314       2.441883789092       1.412518613009
            11    4     -0.279924453037      10.589287680482       1.499991095925
            12    5      5.071900783065       1.574698060931       2.206250402429
            13    6      5.363296139400      11.476422044003       2.265861279676
            14    1      1.488250741855       8.327116753407       6.377715921031
            15    1      0.526080049483       6.542273724088       1.852472804859
            16    1      2.355032633581       8.412777569528       1.428717220956
            17    1      2.261760529027       4.611043770759       1.542107476588
            18    1      4.162379495367      10.631635094440       1.639299124272
            19    1      3.856398432415       2.484246760626       1.739082549228
            20    1      4.911944406028       8.103046095789       2.620934765882
            21    1      4.677388036393       4.802437380894       2.657867592002
            22    2      6.516955192919       6.379073141726       0.885739882561
            23    3      4.001658265295       8.906566107098       1.412518613009
            24    4      3.777358832354       4.124605362477       1.499991095925
            25    5      1.014617497674       8.039380378936       2.206250402429
            26    6      0.847279424791       5.011739725998       2.265861279676
            27    1      8.903405515045      11.112206810970       0.101727257909
            28    1      1.292276207417      12.897049840289       4.626970374082
            29    1     -0.536676376682      11.026545994849       5.050725957984
            30    1     -0.902137701345       1.898915157608       4.937335702353
            31    1      6.229276761533       8.807688469937       4.840144054669
            32    1      6.076524395267       4.025712167740       4.740360629713
            33    1      5.479711850872      11.336277468589       3.858508413059
            34    1      5.255534791289       1.707521547473       3.821575586939
            35    2      3.415967634763       0.130885786641       5.593703296380
            36    3      6.389997991605      10.532757457279       5.066924565932
            37    4      6.155563995328       2.385353565889       4.979452083016
            38    5      0.803738759226      11.399943185441       4.273192776511
            39    6      0.512343402890       1.498219202369       4.213581899265
            40    1      4.387388800436       4.647524492965       0.101727257909
            41    1      5.349559492808       6.432367522284       4.626970374082
            42    1      3.520606908709       4.561863676844       5.050725957984
            43    1      3.613879013264       8.363597475613       4.937335702353
            44    1      1.713260046924       2.343006151931       4.840144054669
            45    1      2.019241109876      10.490394485746       4.740360629713
            46    1      0.963695136263       4.871595150583       3.858508413059
            47    1      1.198251505898       8.172203865478       3.821575586939
            48    2     -0.641315650628       6.595568104646       5.593703296380
            49    3      1.873981276996       4.068075139274       5.066924565932
            50    4      2.098280709937       8.850035883895       4.979452083016
            51    5      4.861022044617       4.935260867436       4.273192776511
            52    6      5.028360117499       7.962901520374       4.213581899265
    """)

    Microcline_Dummy_Elements = ase.io.lammpsdata.read_lammps_data(Microcline_with_Al_Order_Logic_Embedded_File, style="atomic", sort_by_id=True)


    config = ase.build.make_supercell(Microcline_Dummy_Elements, [[x_dim, 0, 0], [0, y_dim, 0], [0, 0, z_dim]])

    return config

# %%
def allUnique(x):
    seen = set()
    return not any(i in seen or seen.add(i) for i in x)


def parse_command_line(args=None):
    parser = argparse.ArgumentParser(description="Process some stuff.")
    parser.add_argument("--Na_fraction", type=float, help="")
    parser.add_argument("--T1O_frac", type=float, nargs='?', default=1.0)
    parser.add_argument("--T1M_frac", type=float, nargs='?', default=0.0)
    parser.add_argument("--T2O_frac", type=float, nargs='?', default=0.0)
    parser.add_argument("--T2M_frac", type=float, nargs='?', default=0.0)

    parser.add_argument('--supercell_dim', type=int, nargs=3, help="")

    parser.add_argument("--sqsGen_NaK", action="store_true")

    args = parser.parse_args()

    return args

# %%

def main():

    args = parse_command_line()

    #T1M_bool = args.T1M_bool
    #T2O_bool = args.T2O_bool
    #T2M_bool = args.T2M_bool
    #T1O_frac = args.T1O_frac

    T1O_frac = args.T1O_frac
    T1M_frac = args.T1M_frac
    T2O_frac = args.T2O_frac
    T2M_frac = args.T2M_frac

    if not (T1O_frac + T1M_frac + T2O_frac + T2M_frac) == float(1):
        raise ValueError("Fractions do not sum to 1")

    T1O_bool = (T1O_frac > 0)
    T1M_bool = (T1M_frac > 0)
    T2O_bool = (T2O_frac > 0)
    T2M_bool = (T2M_frac > 0)

    print("T1O_bool = ", T1O_bool)
    print("T1M_bool = ", T1M_bool)
    print("T2O_bool = ", T2O_bool)
    print("T2M_bool = ", T2M_bool)

    #if T1M_frac > 0:
    #    T1M_bool = True
    #if T2O_frac > 0:
    #    T2O_bool = True
    #if T2M_frac > 0:
    #    T2M_bool = True

    sqsgenerator_NaK = args.sqsGen_NaK

    x_dim, y_dim, z_dim = args.supercell_dim
    Na_fraction = args.Na_fraction
    NaK_SQS_bool = False

    retry = False

    config = build_reference_system(x_dim, y_dim, z_dim)
    print("Built atoms template:", config)

    T1O_site_IDS = []
    T1M_site_IDS = []
    T2O_site_IDS = []
    T2M_site_IDS = []

    elements_initialize = []

    for i in range(len(config.get_chemical_symbols())):
        element = config.get_chemical_symbols()[i]
        
        if element == 'Li':         # Type 3 becomes Li, flag for T1O
            T1O_site_IDS.append(i)
        elif element == 'Be':       # Type 4 -> T1M
            T1M_site_IDS.append(i)
        elif element == 'B':        # Type 5 -> T2O
            T2O_site_IDS.append(i)
        elif element == 'C':        # Type 6 -> T2M
            T2M_site_IDS.append(i)

        # We want to convert to start with ordered Microcline
        if element == 'H':
            elements_initialize.append('O')
        elif element == 'He':
            elements_initialize.append('K')
        elif element == 'Li':
            elements_initialize.append('Al')
        else:
            elements_initialize.append('Si')

    Possible_Al_sites = []
    #Possible_Al_sites.extend(T1O_site_IDS)

    if T1O_bool:
        Possible_Al_sites.extend(T1O_site_IDS)
    if T1M_bool:
        Possible_Al_sites.extend(T1M_site_IDS)
    if T2O_bool:
        Possible_Al_sites.extend(T2O_site_IDS)
    if T2M_bool:
        Possible_Al_sites.extend(T2M_site_IDS)


    config.set_chemical_symbols(elements_initialize)

    element_list = config.get_chemical_symbols()

    ### We first construct the Na/K distribution
    if sqsgenerator_NaK:
        print("""sqsGenerator is a software developed by Dominik Gehringer
                 and is available at https://github.com/dgehringer/sqsgenerator
                 see also "Models of configurationally-complex alloys made simple": https://doi.org/10.1016/j.cpc.2023.108664""")

        N_Na = int(element_list.count('K')*Na_fraction)
        N_K = element_list.count('K') - N_Na

        sqs_config = from_ase_atoms(config)  # Convert to sqsgenerator structure

        try:
            configuration = dict(
                structure=sqs_config,
                iterations=1e7, # 1e9 is in guide
                shell_weights={1: 1.0,
                            2: 1/2,
                            3: 1/3,
                            4: 1/4,
                            5: 1/5,
                            6: 1/6,
                            7: 1/7,
                            8: 1/8,
                            9: 1/9,
                            10: 1/10
                            },
                which='K',
                composition=dict(Na=N_Na, K=N_K)
            )

            results = sqsgenerator.public.sqs_optimize(configuration, pass_structure=True, structure_format='ase', make_structures=True)

        except sqsgenerator.settings.exceptions.BadSettings:
            # system probably too small for shell 10
            configuration = dict(
                structure=sqs_config,
                iterations=1e7, # 1e9 is in guide
                which='K',
                composition=dict(Na=N_Na, K=N_K)
            )

            results = sqsgenerator.public.sqs_optimize(configuration, pass_structure=True, structure_format='ase', make_structures=True)

        print("sqsGenerator results for Na/K distribution:")
        print(results)
        config_return = results[0][next(iter(results[0]))]['structure']
        print(config_return)

        config = config_return
        #result = sqsgenerator.public.to_ase_atoms(results)
    else:
        ### Just change Na_fraction * K atoms to Na atoms randomly
        K_to_change = random.sample(range(element_list.count('K')), int(element_list.count('K')*Na_fraction))

        element_transform = []
        K_Count = -1
        for i in range(len(element_list)):
            if element_list[i] == 'K':
                K_Count += 1
                if K_Count in K_to_change:
                    element_transform.append('Na')
                else:
                    element_transform.append('K')
            else:
                element_transform.append(element_list[i])

        config.set_chemical_symbols(element_transform)

    element_list = config.get_chemical_symbols()

    ########################################################
    ###### Al and Si disorder is done in three stages ######
    ########################################################
    ###### 1. Random distribution of Al and Si atoms
    ###### 2. Random distribution of Al and Si atoms + Löwenstein rule is fulfilled
    ###### 3. Random distribution of Al and Si atoms + Löwenstein rule is fulfilled + T fractions are fulfilled

    ### The algorithm has no problems with the concentrations that one might see in Nature like:
            # Any T1O + T1M                     (such as --T1O_frac 0.5 --T1M_frac 0.5 ~~ Orthoclase (kinda))
            # Any T1O + T1M + T2O + T2M         (such as --T1O_frac 0.25 --T1M_frac 0.25 --T2O_frac 0.25 --T2M_frac 0.25 ~~ Sanidine)
            # Just T1O                          (such as --T1O_frac 1.0 ~~ Microcline)
    ### For strange concentrations like (T1O + T2M) or (T2O + T1M) the algorithm has to take some roundabout ways
    ### (basically optimizer steps in the wrong direction to get over a hill) that might take significant longer 
    ### to converge since its not tweaked for these cases.



    ###### 1. Random distribution of Al and Si atoms ######

    Al_to_change = random.sample(Possible_Al_sites, int(element_list.count('Al')))

    print("allUnique() Sanity Checks passed?:", allUnique(Possible_Al_sites), allUnique(Al_to_change))

    element_transform = []

    for i in range(len(element_list)):
        if i in Al_to_change and not(i in Possible_Al_sites):
            raise ValueError("Al_to_change" + str(i) + " not in Possible_Al_sites")
        if element_list[i] == 'Al':
            element_transform.append('Si')
        else:
            element_transform.append(element_list[i])

    for id in Al_to_change:
        element_transform[id] = 'Al'

    config.set_chemical_symbols(element_transform)


    ###### 2. Random distribution of Al and Si atoms + Löwenstein rule is fulfilled ######

    cutoff = [1.7]*len(element_list)
    neighbours = ase.neighborlist.NeighborList(cutoff, skin=0, self_interaction=False, bothways=True)
    neighbours.update(config)

    for k in range(2000):
        AlNr = 0
        SiNr = 0
        Possible_Si_sites_w_4_Si = []
        Al_sites_w_Al_n = []
        for i in range(len(config.get_chemical_symbols())):
            if config.get_chemical_symbols()[i] == "Al":
                AlNr += 1
                no_Al = True
                Si_Al_sanity_int = 0
                for n in neighbours.get_neighbors(i)[0]:
                    if config.get_chemical_symbols()[n] == 'Al':
                        no_Al = False
                        Si_Al_sanity_int += 1
                    elif config.get_chemical_symbols()[n] == 'Si':
                        Si_Al_sanity_int += 1
                if not no_Al:
                    Al_sites_w_Al_n.append(i)

                if Si_Al_sanity_int != 4:
                    raise ValueError("Al " + str(i) + " has " + str(Si_Al_sanity_int) + " neighbours.")

            if config.get_chemical_symbols()[i] == "Si" and i in Possible_Al_sites:
                SiNr += 1
                no_Al = True
                Si_Al_sanity_int = 0

                for n in neighbours.get_neighbors(i)[0]:
                    if config.get_chemical_symbols()[n] == 'Al':
                        no_Al = False
                        Si_Al_sanity_int += 1
                    elif config.get_chemical_symbols()[n] == 'Si':
                        Si_Al_sanity_int += 1

                if no_Al:
                    Possible_Si_sites_w_4_Si.append(i)
                if Si_Al_sanity_int != 4:
                    raise ValueError("Si " + str(i) + " has " + str(Si_Al_sanity_int) + " neighbours.")

        element_list = config.get_chemical_symbols()
        #print("Al_sites_w_Al_n:", Al_sites_w_Al_n)
        #print("Possible_Si_sites_w_4_Si", Possible_Si_sites_w_4_Si)

        Al_sites_w_Al_n = random.sample(Al_sites_w_Al_n, int(len(Al_sites_w_Al_n)/2))
        for Al in Al_sites_w_Al_n:
            if not Possible_Si_sites_w_4_Si:
                print("No Si sites available for Al sites. This should only happen for strange T concentrations (T1O + T2M) or (T2O + T1M)")
                print("Allow temporary all sites")
                Possible_Al_sites = []
                Possible_Al_sites.extend(T1O_site_IDS)
                Possible_Al_sites.extend(T1M_site_IDS)
                Possible_Al_sites.extend(T2O_site_IDS)
                Possible_Al_sites.extend(T2M_site_IDS)
                T1O_bool = True
                T1M_bool = True
                T2O_bool = True
                T2M_bool = True

                break

            local_Si = random.randrange(len(Possible_Si_sites_w_4_Si))
            Si = Possible_Si_sites_w_4_Si.pop(local_Si)
            element_list[Al] = 'Si'
            element_list[Si] = 'Al'
        config.set_chemical_symbols(element_list)
        if len(Al_sites_w_Al_n) == 0:
            break
    if not len(Al_sites_w_Al_n) == 0:
        print("Al Sites not converged. Retry please")
        exit()


    ###### 3. Random distribution of Al and Si atoms + Löwenstein rule is fulfilled + T fractions are fulfilled ######

    correct_fraction_bool = False
    iteration = 0

    while not correct_fraction_bool:
        element_list = config.get_chemical_symbols()

        iteration += 1
        print("##### T FRAC OPT ITERATION:", iteration)

        T1O_fraction = 0
        T1M_fraction = 0
        T2O_fraction = 0
        T2M_fraction = 0

        Possible_Si_T1O_sites_w_4_Si = []
        Possible_Si_T1M_sites_w_4_Si = []
        Possible_Si_T2O_sites_w_4_Si = []
        Possible_Si_T2M_sites_w_4_Si = []

        Al_T1O_sites = []
        Al_T1M_sites = []
        Al_T2O_sites = []
        Al_T2M_sites = []

        for i in range(len(element_list)):
            if element_list[i] == 'Al':
                if i in T1O_site_IDS:
                    T1O_fraction += 1
                    Al_T1O_sites.append(i)
                elif i in T1M_site_IDS:
                    T1M_fraction += 1
                    Al_T1M_sites.append(i)
                elif i in T2O_site_IDS:
                    T2O_fraction += 1
                    Al_T2O_sites.append(i)
                elif i in T2M_site_IDS:
                    T2M_fraction += 1
                    Al_T2M_sites.append(i)

            elif element_list[i] == 'Si':
                no_Al = True
                for n in neighbours.get_neighbors(i)[0]:
                    if config.get_chemical_symbols()[n] == 'Al':
                        no_Al = False
                if no_Al:
                    if i in T1O_site_IDS:
                        Possible_Si_T1O_sites_w_4_Si.append(i)
                    elif i in T1M_site_IDS:
                        Possible_Si_T1M_sites_w_4_Si.append(i)
                    elif i in T2O_site_IDS:
                        Possible_Si_T2O_sites_w_4_Si.append(i)
                    elif i in T2M_site_IDS:
                        Possible_Si_T2M_sites_w_4_Si.append(i)
        for i in Possible_Si_T1O_sites_w_4_Si:
            for n in neighbours.get_neighbors(i)[0]:
                if config.get_chemical_symbols()[n] == 'Al':
                    raise ValueError("Possible_Si " + str(i) + " has Al neighbour.")
        for i in Possible_Si_T1M_sites_w_4_Si:
            for n in neighbours.get_neighbors(i)[0]:
                if config.get_chemical_symbols()[n] == 'Al':
                    raise ValueError("Possible_Si " + str(i) + " has Al neighbour.")
        for i in Possible_Si_T2O_sites_w_4_Si:
            for n in neighbours.get_neighbors(i)[0]:
                if config.get_chemical_symbols()[n] == 'Al':
                    raise ValueError("Possible_Si " + str(i) + " has Al neighbour.")
        for i in Possible_Si_T2M_sites_w_4_Si:
            for n in neighbours.get_neighbors(i)[0]:
                if config.get_chemical_symbols()[n] == 'Al':
                    raise ValueError("Possible_Si " + str(i) + " has Al neighbour.")

        T1O_fraction=T1O_fraction/element_list.count('Al')
        T1M_fraction=T1M_fraction/element_list.count('Al')
        T2O_fraction=T2O_fraction/element_list.count('Al')
        T2M_fraction=T2M_fraction/element_list.count('Al')

        print("Al#  = ", element_list.count('Al'))
        print("T1O_fraction = ", T1O_fraction)
        print("T1M_fraction = ", T1M_fraction)
        print("T2O_fraction = ", T2O_fraction)
        print("T2M_fraction = ", T2M_fraction)

        if (T1O_bool + T1M_bool + T2O_bool + T2M_bool) == 2:
            #### If only a single setting is present we remove one Al before moving Al since its too crowded.
            if (T1O_fraction > T1O_frac or T1M_fraction > T1M_frac or T2O_fraction > T2O_frac or T2M_fraction > T2M_frac):
                Al_sites = [Al_T1O_sites, Al_T1M_sites, Al_T2O_sites, Al_T2M_sites]
                T_site_ids = [T1O_site_IDS, T1M_site_IDS, T2O_site_IDS, T2M_site_IDS]
                Al_atom_triage = [T1O_fraction - T1O_frac, T1M_fraction - T1M_frac, T2O_fraction - T2O_frac, T2M_fraction - T2M_frac]
                index_too_many_Als = Al_atom_triage.index(max(Al_atom_triage))
                index_too_few_Als = Al_atom_triage.index(min(Al_atom_triage))

                Al_possible_change = Al_sites[index_too_many_Als]

                # Already remove one Al
                Al_to_change = random.randrange(len(Al_possible_change))
                Al = Al_possible_change[Al_to_change]
                element_list[Al] = 'Si'
                config.set_chemical_symbols(element_list)

                Possible_other_Si_sites_w_4_Si = []

                Al_T1O_sites = []
                Al_T1M_sites = []
                Al_T2O_sites = []
                Al_T2M_sites = []

                Changed_sites = [Al, *neighbours.get_neighbors(Al)[0]]
                print("Changed sites:", Changed_sites)
                for i in Changed_sites:
                    if element_list[i] == 'Si':
                        no_Al = True
                        for n in neighbours.get_neighbors(i)[0]:
                            if config.get_chemical_symbols()[n] == 'Al':
                                no_Al = False
                        if no_Al:
                            if i in T_site_ids[index_too_many_Als]:
                                Possible_Si_T1O_sites_w_4_Si.append(i)
                            elif i in T_site_ids[index_too_few_Als]:
                                Possible_other_Si_sites_w_4_Si.append(i)

                if Possible_other_Si_sites_w_4_Si:
                    Si_possible_change = Possible_other_Si_sites_w_4_Si
                else:
                    Si_possible_change = Possible_Si_T1O_sites_w_4_Si

                Si_to_change = random.randrange(len(Si_possible_change))
                Si = Si_possible_change[Si_to_change]
                element_list[Si] = 'Al'

                config.set_chemical_symbols(element_list)
            else:
                correct_fraction_bool = True

        else:
            if (T1O_fraction > T1O_frac or T1M_fraction > T1M_frac or T2O_fraction > T2O_frac or T2M_fraction > T2M_frac):

                Al_sites = [Al_T1O_sites, Al_T1M_sites, Al_T2O_sites, Al_T2M_sites]
                Al_atom_triage = [T1O_fraction - T1O_frac, T1M_fraction - T1M_frac, T2O_fraction - T2O_frac, T2M_fraction - T2M_frac]
                
                ## Maybe if we choose the index to move/remove randomly, weighted by the triage we would get good convergence 
                ## also for strange concentrations
                index_too_many_Als = Al_atom_triage.index(max(Al_atom_triage))
                index_too_few_Als = Al_atom_triage.index(min(Al_atom_triage))

                Al_possible_change = Al_sites[index_too_many_Als]
                Si_sites = [Possible_Si_T1O_sites_w_4_Si, Possible_Si_T1M_sites_w_4_Si, Possible_Si_T2O_sites_w_4_Si, Possible_Si_T2M_sites_w_4_Si]

                minmax_bool = True
                if minmax_bool:
                    Si_possible_change = Si_sites[index_too_few_Als]
                
                while not Si_possible_change:
                    ####### This should only happen for a strange concentrations.
                    ####### The algorithm should still converge but might take much longer

                    print("No free Sis in", ["T1O", "T1M", "T2O", "T2M"][index_too_few_Als], "Lateral move")
                    Si_possible_change = random.choice(Si_sites)

                Al_to_change = random.randrange(len(Al_possible_change))
                Al = Al_possible_change[Al_to_change]
                Si_to_change = random.randrange(len(Si_possible_change))
                Si = Si_possible_change[Si_to_change]

                print(Al_possible_change, Al)
                print(Si_possible_change, Si)

                element_list[Al] = 'Si'
                element_list[Si] = 'Al'
                config.set_chemical_symbols(element_list)

            else:
                correct_fraction_bool = True

    for i in range(len(config.get_chemical_symbols())):
        if config.get_chemical_symbols()[i] == "Al":
            AlNr += 1
            no_Al = True
            for n in neighbours.get_neighbors(i)[0]:
                if config.get_chemical_symbols()[n] == 'Al':
                    no_Al = False
                    print("ERROR: Al " + str(Al) +  " has Al neighbours. EXIT")
                    exit()

    print("Al#  = ", element_list.count('Al'))
    print("T1O_fraction = ", T1O_fraction)
    print("T1M_fraction = ", T1M_fraction)
    print("T2O_fraction = ", T2O_fraction)
    print("T2M_fraction = ", T2M_fraction)


    ase.io.lammpsdata.write_lammps_data('ASE_2x2x3_microcline_pristine.lmp', config, specorder=['O', 'K', 'Na', 'Al', 'Si'])


if __name__ == "__main__":
    main()


# %%