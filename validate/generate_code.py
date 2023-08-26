import sys
import os

datasets = [
    'Agrafiotis_ABCD',
    #'BindingDB_exact',
    #'BindingDB_similarity',
    #'BindingDB_substructure',
    'Hicks_and_Jochum',
    'RDKit_smarts',
    'RDMACCS',
    'Rarey_smarts',
]

part_size = 50


def exclude_smarts(smarts):
    if '.' in smarts:
        return True, 'components'
    if '$(' in smarts:
        return True, 'recursive'
    for s in ['/', '\\', '@?', '@@', '@]', '@H']:
        if s in smarts:
          return True, 'stereo'

    parse_errors = [
    ]

    if smarts in parse_errors:
        return True, 'parse error'


    return False, None



def read_smarts(sqc_path, dataset):
    with open(os.path.join(sqc_path, dataset, '{}.dat'.format(dataset)), 'r') as f:
        return [line.strip().split()[0] for line in f.readlines()]


def generate_dataset_files(dataset: str, smarts: list, part: int = 0):
    offset = part * part_size
    part_name = '{}_part_{}'.format(dataset, part + 1)

    with open('{}.cpp'.format(part_name), 'w') as f:
        f.write('#include "Validate.hpp"\n')
        f.write('\n')
        f.write('void {}_part_{}(OpenBabel::OBMol &mol)\n'.format(dataset, part + 1))
        f.write('{\n')
        stop = min(offset + part_size, len(smarts))
        f.write('    // SMARTS {} - {}\n'.format(offset + 1, stop))
        for i in range(offset, stop):
            exclude, reason = exclude_smarts(smarts[i])
            comment = '//' if exclude else ''
            reason = ' // FIXME: {}'.format(reason) if exclude else ''
            f.write('    {}validate_contains<"{}">(mol);{}\n'.format(comment, smarts[i], reason))
        f.write('}\n')

    yield part_name

    if offset + part_size < len(smarts):
        yield from generate_dataset_files(dataset, smarts, part + 1)


def generate_part_files(sqc_path):
    for dataset in datasets:
        smarts = read_smarts(sqc_path, dataset)
        print('{}: {}'.format(dataset, len(smarts)))
        yield from generate_dataset_files(dataset, smarts)


def generate_main(part_names):
    with open('ValidateCTSmarts.cpp', 'w') as f:
        f.write('#include <Kitimar/OpenBabel/OpenBabel.hpp>\n')
        f.write('#include "../test/TestData.hpp"\n')
        f.write('#include <gtest/gtest.h>\n')
        f.write('\n')
        f.write('using namespace Kitimar;\n')
        f.write('\n')
        for part_name in part_names:
            f.write('void {}(OpenBabel::OBMol &mol);\n'.format(part_name))
        f.write('\n')
        f.write('TEST(ValidateCTSmarts, Sqc)\n')
        f.write('{\n')
        f.write('    OpenBabelSmilesMolSource source{chembl_smi_filename("100K")};\n')
        f.write('    auto i = 0;\n')
        f.write('    for (auto mol : source.molecules()) {\n')
        f.write("        std::cout << "Molecule: " << ++i << " -- " << writeSmiles(mol) << '\\n';\n")
        for part_name in part_names:
            f.write('        {}(mol);\n'.format(part_name))
        f.write('    }\n')
        f.write('}\n')


def generate_cmake(part_names):
    with open('ValidateCTSmarts.cmake', 'w') as f:
        f.write('set(ValidateCTSmarts_PART_SRCS\n')
        for part_name in part_names:
            f.write('    {}.cpp\n'.format(part_name))
        f.write(')\n')


def generate_code(sqc_path):
    part_names = list(generate_part_files(sqc_path))
    generate_main(part_names)
    generate_cmake(part_names)
    
        
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: {} <path to sqc>'.format(sys.argv[0]))
        raise SystemExit(1)

    generate_code(sys.argv[1])
