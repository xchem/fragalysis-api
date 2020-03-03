import luigi
from fragalysis_api import Align, set_up
import os
import glob
import pandas as pd


def outlist_from_align(input_dir, output_dir):
    outlist = []
    for f in glob.glob(os.path.join(input_dir, "*.pdb")):
        name = f.split('/')[-1].replace('.pdb', '')
        out = os.path.join(output_dir, 'tmp', f'{name}_bound.pdb')
        outlist.append(out)
    return outlist


class AlignTarget(luigi.Task):
    input_dir = luigi.Parameter()
    output_dir = luigi.Parameter()

    def requires(self):
        pass

    def output(self):
        outlist = outlist_from_align(input_dir=self.input_dir, output_dir=self.output_dir)
        return [luigi.LocalTarget(o) for o in outlist]

    def run(self):
        structure = Align(self.input_dir, pdb_ref="")
        structure.align(os.path.join(self.output_dir, "tmp"))


class ProcessAlignedPDB(luigi.Task):
    input_dir = luigi.Parameter()
    input_file = luigi.Parameter()
    target_name = luigi.Parameter()
    output_dir = luigi.Parameter()

    def requires(self):
        return AlignTarget(input_dir=self.input_dir, output_dir=self.output_dir)

    def output(self):
        return luigi.LocalTarget(os.path.join(self.output_dir, 'mol_dict.csv'))

    def run(self):
        new = set_up(
            target_name=self.target_name,
            infile=self.input_file,
            out_dir=self.output_dir,
        )

        df = pd.DataFrame.from_dict(new.mol_dict)
        df.to_csv(self.output().path)


class BatchProcessAlignedPDB(luigi.Task):
    input_dir = luigi.Parameter()
    target_name = luigi.Parameter()
    output_dir = luigi.Parameter()

    def requires(self):
        self.aligned_list = outlist_from_align(self.input_dir, self.output_dir)
        return [
            ProcessAlignedPDB(
                target_name=self.target_name, input_file=i, output_dir=self.output_dir,
                input_dir=self.input_dir
            )
            for i in self.aligned_list
        ]

    def output(self):
        return luigi.LocalTarget(os.path.join(self.output_dir, 'aligned_list.pdb'))

    def run(self):
        with open(self.output().path, 'w') as f:
            f.write('\n'.join(self.aligned_list))
        f.close()



class BatchConvertAligned(luigi.Task):
    search_directory = luigi.Parameter()
    output_directory = luigi.Parameter()

    def requires(self):
        in_lst = [os.path.abspath(f.path) for f in os.scandir(self.search_directory) if f.is_dir()]
        out_lst = []
        target_names = []

        for f in in_lst:
            out = os.path.join(os.path.abspath(self.output_directory), f.split('/')[-1])
            out_lst.append(out)
            target_names.append(f.split('/')[-1]).upper()

        return [BatchProcessAlignedPDB(input_dir=i, output_dir=self.output_directory, target_name=t)
                for (i, t) in list(zip(in_lst, target_names))]

    def output(self):
        return luigi.LocalTarget(os.path.join(self.search_directory, 'dir_list.txt'))

    def run(self):
        lst = [os.path.abspath(f.path) for f in os.scandir(self.search_directory) if f.is_dir()]
        lst_str = '\n'.join([f for f in lst])
        with open(self.output().path, 'w') as w:
            w.write(lst_str)
        w.close()