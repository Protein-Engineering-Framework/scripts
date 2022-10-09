# Niklas Siedhoff, 16 Aug 2022, drawing residue distances
# Requires PyMOL version 2.x (or higher) 
# Run me in PyMOL command-line after opening 2v1a.pse with 
#       run draw_dash_lines_PyMOL.py

from pymol.cmd import distance, select, color, show, label, alter

coupled_pairs = {
    # comment/outcomment variant-fitness pairs as desired for plotting
    # Zayner 2012
    #'W491Y/L408W': -0.0163904161881694,
    #'T406A/T407A': -0.160296991872449,
    #'E412H/Q436H': -0.317420411852151,
    # Zayner 2013
    #'D432A/Q436A': -0.0397671268714877,
    'S441A/E443A': 0.0263289387223491,
    #'D432A/Q436A/E443A': -0.0222763947111523,
    #'D432A/Q436A/T406A/T407A': -0.11069829749369,
    #'S441A/E443A/T406A/T407A': -0.0579919469776868,
    # Zayner 2014
    #'N414A/Q513H': -1.60745502321467,
    'N414L/Q513A': 1.37026858207418,
    'N414A/Q513A': 1.40978706133298,
    #'C450V/F494C': -0.0696359281413944, ## makes no sense that a C450 mutant shows a recovery, delete
    #'C450V/Q513C': 44, ## makes no sense that a C450 mutant shows a recovery, delete
    #'N414G/C450V/Q513C': 13, ## makes no sense that a C450 mutant shows a recovery, delete
    # Zoltowski 2009
    'V416I/L496I': 1.09540614735826,
    # Christie 2007
    #'K413R/I427V': -1.07918124604762,
    #'K413R/I427L': -0.499397649430815,
    # Hemmer 2022 Dataset, 1st round
    #'I427T/E475T': -1.31175386105575,
    'Q513R/G528K': 0.179607832778518,
    'Q513R/G528R': 0.3556990918342,
    'K413C/N414C': 0.62526224640906,
    #'L446M/D501W': -0.135662602000073,
    #'L446Q/D501Y': -0.466655821041497,
    #'L446S/D501Y': -0.767685816705479,
    # Hemmer 2022 Dataset, 2nd round
    #'E475T/G528E': -0.10763387839983,
    'N414G/G528E': 0.772822416878577,
    'R410P/G528A': 0.165367393663908,
    'R410P/G528E': 0.150644136843202,
    'R410P/G528R': 0.030668819766452,
    'K413A/N414G': 0.810462017217072,
    'N414G/L446E': 0.0404286570556082,
    'N414G/L446S': 1.10321948691506,
    'R442L/D515L': 0.172545978291032,
    #'K413A/D501G': -0.0217192496932363,
    'D501G/G528E': 0.0499739749618386,
    'N414A/V416A': 0.010465433678165,
    'N414A/V416L': 3.16365580035721,
    'N414L/V416L': 3.40013768882482,
    'N414L/L514A': 1.40089584057146,
    'N414L/H495L': 1.06207728401808,
    'N414L/Q513L': 2.04110370135124,
    # Hemmer 2022 Dataset, 3rd round
    #'T438V/E475T': -0.382334935341462,
    'V463W/E475T': 0.0499739749618386,
    #'N414D/V416T': -1.01072386539177,
    #'I427T/L446M/E475T': -2.01072386539177,
    #'N414D/F415E/V416T': -0.094269916841848
}

def color_fitness(fitness):
    if fitness < -0.08:  # lower threshold for fast recovery, red (change here)
        return 'red'
    elif fitness > 0.10:  # upper threshold for slow recovery, blue (change here)
        return 'blue'
    else:  # else similar to WT, orange
        return 'orange'


target = '2v1a'
print_distance = 0  # 0: do not show distance, 1: show distance
print_label = 0  # 0: do not show residue name, 1: show residue name
verbose = 1  # 0: do not print commands, 1: print commands
vdW_radius_size = 1.5  # change the size of the vdW ball at CA atoms
# counter for dist command, here in this PyMOL scene starting with 'dist01' --> counter+=0 below
# else, e.g., add 3 (counter = 3) if starting with dist04
counter = 0
paired_pos = []
if verbose:
    print('#'*60 + f'\nInput: Totally {len(coupled_pairs)} variant-dark recovery pairs.')
    print('Executing commands for each variant-fitness pair...\n' + '#'*60 + '\n')
for i, pair in enumerate(coupled_pairs.items()):
    pair_variant = pair[0]
    recovery_time = pair[1]
    if verbose:
        print(f'######### {i+1}: Variant {pair_variant}, Dark rec.: {recovery_time} s #########')
    color_ = color_fitness(recovery_time)
    pair_split = pair_variant.split('/')
    positions = []
    for single in pair_split:
        position = single[1:-1]
        positions.append(position)
        show('spheres', f'(/{target}///{position}/CA)')
        alter(f'(/{target}///{position}/CA)', f'vdw="{vdW_radius_size}"')
        if verbose:
            print(f'show spheres, (/{target}///{position}/CA)')
            print(f'alter (/{target}///{position}/CA), vdw={vdW_radius_size}')
        #color(f'{color_}', f'resi {position}')  # just color positions white
        #if verbose:
        #    print(f'show sticks, resi {position}')
        #    print(f'color {color_}, resi {position}')
        if print_label:
            label(f'(/{target}///{position}/CA)', f'{position}')
            if verbose:
                print(f'label (/{target}///{position}/CA), {position}')
    distance_cmds = []
    for j, pos in enumerate(positions):
        if j < len(positions)-1:
            counter += 1
            paired_pos.append((positions[j], positions[j + 1]))
            if paired_pos.count((positions[j], positions[j + 1])) == 1:
                distance(f'(/{target}///{positions[j]}/CA)', f'(/{target}///{positions[j+1]}/CA)', label=print_distance)
                if verbose:
                    print(f'distance (/{target}///{positions[j]}/CA), (/{target}///{positions[j + 1]}/CA)')
            elif paired_pos.count((positions[j], positions[j + 1])) == 2:
                distance(f'(/{target}///{positions[j]}/C)', f'(/{target}///{positions[j+1]}/C)', label=print_distance)
                if verbose:
                    print(f'distance (/{target}///{positions[j]}/C), (/{target}///{positions[j + 1]}/C)')
            else:
                distance(f'(/{target}///{positions[j]}/N)', f'(/{target}///{positions[j + 1]}/N)', label=print_distance)
                print(f'distance (/{target}///{positions[j]}/N), (/{target}///{positions[j+1]}/N)')
            color(f'{color_}', f'dist{counter:02d}')
            if verbose:
                print(f'color {color_}, dist{counter:02d}')
if verbose:
    print('\n' + '#'*60 + '\nDone!\n' + '#'*60)

