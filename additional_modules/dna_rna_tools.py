def valid_sequence(seqs):
    for nucl_acid in seqs:
        t_in_nucl_acid = "T" in nucl_acid or "t" in nucl_acid
        u_in_nucl_acid = "U" in nucl_acid or "u" in nucl_acid
        if t_in_nucl_acid and u_in_nucl_acid:
            return False
    return True


def seq_transcr(seq):
    return seq.replace("T", "U").replace("t", "u")


def seq_reverse(seq):
    return seq[::-1]


def seq_compl(seq):
    dna_complement_dict = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "a": "t",
        "t": "a",
        "c": "g",
        "g": "c",
        "U": "A",
        "u": "a",
    }

    return "".join(
        dna_complement_dict[base] if base in dna_complement_dict else base
        for base in seq
    )


def seq_rev_compl(seq):
    return seq_reverse(seq_compl(seq))


if __name__ == "__main__":
    seqs = input("Введите последовательность нуклеиновой кислоты: ").split()
    procedure = input("Введите тип операции: ")
    result = run_dna_rna_tools(*seqs, procedure)
    print(result)
