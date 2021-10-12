from dataclasses import dataclass
from collections import defaultdict
from itertools import combinations, product
from enum import Enum


class Sex(Enum):
    UNKNOWN = 0
    MALE = 1
    FEMALE = 2


class Phenotype(Enum):
    MISSING_ONE = -9
    MISSING_TWO = 0
    CONTROL = 1
    CASE = 2


class Pedigree:
    def __init__(self, fh):
        self.families = {}
        for line in fh:
            (
                family_id,
                sample_id,
                father_id,
                mother_id,
                sex,
                phenotype,
            ) = line.strip().split("\t")
            if family_id not in self.families:
                self.families[family_id] = Family(family_id)
            self.families[family_id].add_sample(
                sample_id, father_id, mother_id, sex, phenotype
            )
        self.validate_ped()

    def get_num_samples(self):
        nsamples = 0
        for family in self.families.values():
            nsamples += len(family.get_sample_ids())
        return nsamples

    def get_all_sample_ids(self):
        sample_ids = set()
        for family in self.families.values():
            sample_ids |= family.get_sample_ids()
        return sample_ids

    def get_num_families(self):
        return len(self.families)

    def get_intra_family_relationships(self):
        child_parent_pairs, sibling_pairs, others = set(), set(), set()
        for family in self.families.values():
            cp_pairs, s_pairs, o_pairs = family.get_relationship_pairs()
            child_parent_pairs |= cp_pairs
            sibling_pairs |= s_pairs
            others |= o_pairs
        return child_parent_pairs, sibling_pairs, others

    def get_inter_family_relationships(self):
        for family_a, family_b in combinations(sorted(self.families), 2):
            for sample_a, sample_b in product(
                self.families[family_a].get_sample_ids(),
                self.families[family_b].get_sample_ids(),
            ):
                yield (sample_a, sample_b)

    def validate_ped(self):
        samples_so_far = set()
        for family in self.families.values():
            family.validate_family()
            family_samples = family.get_sample_ids()
            overlap = samples_so_far & family_samples
            if overlap:
                raise ValueError(
                    f"Samples {', '.join(overlap)} observed multiple times."
                )
            samples_so_far |= family_samples

    def __str__(self):
        return "\n".join((str(family) for family in self.families.values()))


class Family:
    def __init__(self, family_id):
        self.family_id = family_id
        self.samples = {}

    def add_sample(self, sample_id, father_id, mother_id, sex, phenotype):
        if sample_id in self.samples:
            raise ValueError(
                f"{sample_id} is already present in family {self.family_id}."
            )
        self.samples[sample_id] = Sample(
            self.family_id, sample_id, father_id, mother_id, sex, phenotype
        )

    def get_sample_ids(self):
        return set(self.samples.keys())

    def get_relationship_pairs(self):
        others = set(combinations(sorted(self.samples), 2))
        child_parent_pairs = self.get_child_parent_pairs()
        others -= child_parent_pairs
        sibling_pairs = self.get_sibling_pairs()
        others -= sibling_pairs
        num_samples = len(self.samples)
        expected_num_pairs = int(num_samples * (num_samples - 1) / 2)
        num_pairs = len(child_parent_pairs) + len(sibling_pairs) + len(others)
        if expected_num_pairs != num_pairs:
            raise ValueError(
                f"Error with family {self.family_id}: found "
                f"{num_pairs} pairs but should have "
                f"{expected_num_pairs} for {num_samples} samples."
            )
        return child_parent_pairs, sibling_pairs, others

    def get_child_parent_pairs(self):
        pairs = []
        for sample_id, sample in self.samples.items():
            if sample.father_id:
                pairs.append(tuple(sorted((sample_id, sample.father_id))))
            if sample.mother_id:
                pairs.append(tuple(sorted((sample_id, sample.mother_id))))
        return set(pairs)

    def get_sibling_pairs(self):
        """
        N.B. A PED file doesn't actually have a way to distinguish a sibling
        from a replicate or identical twin.
        """
        samples_by_parents = defaultdict(list)
        for sample_id, sample in self.samples.items():
            if sample.father_id and sample.mother_id:
                samples_by_parents[(sample.father_id, sample.mother_id)].append(
                    sample_id
                )
        pairs = []
        for samples in samples_by_parents.values():
            pairs.extend(combinations(sorted(samples), 2))
        return set(pairs)

    def validate_family(self):
        for sample_id, sample in self.samples.items():
            if sample.father_id:
                if sample.father_id in self.samples:
                    if self.samples[sample.father_id].sex != Sex.MALE:
                        raise ValueError(
                            f"Sample {sample_id} has father {self.samples[sample.father_id]} "
                            "which is not labeled as male."
                        )
                else:
                    raise ValueError(
                        f"Sample {sample_id} has father {sample.father_id} "
                        "which is not present."
                    )
            if sample.mother_id:
                if sample.mother_id in self.samples:
                    if self.samples[sample.mother_id].sex != Sex.FEMALE:
                        raise ValueError(
                            f"Sample {sample_id} has mother {self.samples[sample.mother_id]} "
                            "which is not labeled as female."
                        )
                else:
                    raise ValueError(
                        f"Sample {sample_id} has mother {sample.mother_id} "
                        "which is not present."
                    )

    def __str__(self):
        return "\n".join((str(sample) for sample in self.samples.values()))


@dataclass
class Sample:
    family_id: str
    sample_id: str
    father_id: str
    mother_id: str
    sex: Sex
    phenotype: Phenotype

    def __post_init__(self):
        if self.father_id == "0":
            self.father_id = None
        if self.mother_id == "0":
            self.mother_id = None

        self.sex = Sex(int(self.sex))
        self.phenotype = Phenotype(int(self.phenotype))

    def __str__(self):
        return "\t".join(
            (
                self.family_id,
                self.sample_id,
                self.father_id if self.father_id else "0",
                self.mother_id if self.mother_id else "0",
                str(self.sex.value),
                str(self.phenotype.value),
            )
        )
