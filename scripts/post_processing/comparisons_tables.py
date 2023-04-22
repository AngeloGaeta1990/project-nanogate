import pandas as pd
import argh
# # inputs
# minimap2_stats_filename = "minimap2_summary_table.csv"
# minimap2_combinatorial_filename = "minimap2_combinatorial_summary_table.csv"
# nanogate_positives_stats_filename = "jaccard_summary_table.csv"
# nanogate_time_filename = "nanogate_summary_table.csv"
#
# #outputs
# time_comparison_filename = "time_comparison.csv"
# true_positives_comparison_filename = "true_postives_comparison.csv"


def add_method_row(table, method):
    """
    function to add a row to stats table specifying the method : nanogate, combinatorial, minimap2
    :param table:
    :param method:
    :return:
    """
    row_number = len(table)
    table_method = [method] * row_number
    table["Method"] = table_method
    return table

def comparison_table(

                     # inputs
                     minimap2_stats_filename: 'minimap2 summary stats',
                     nanogate_positives_stats_filename: 'nanogate positives stats filename',
                     nanogate_time_filename: 'nanogate time filename',
                     minimap2_combinatorial_filename: 'minimap2 combinatorial filename',
                     minimap2_constructs_stats_filename: 'minimap2 stats per construct',
                     nanogate_constructs_stats_filename: 'nanogate stats per construct',
                     combinatorial_stats_filename: 'combinatorial stats per construct',
                     # outputs
                     time_comparison_filename: 'time comparison filename',
                     true_positives_comparison_filename: 'true positives comparison filename',
                     constructs_comparison_filename: 'constructs table comparison filename'):
    """
    method to join nanogate and minimap2 tables
    """

    """
    Reading tables
    """
    minimap2_stats = pd.read_csv(minimap2_stats_filename)
    nanogate_positives = pd.read_csv(nanogate_positives_stats_filename)
    nanogate_time = pd.read_csv(nanogate_time_filename)
    combinatorial_stats = pd.read_csv(minimap2_combinatorial_filename)
    minimap2_constructs_table = pd.read_csv(minimap2_constructs_stats_filename)
    nanogate_constructs_table = pd.read_csv(nanogate_constructs_stats_filename)
    combinatorial_constructs_table = pd.read_csv(combinatorial_stats_filename)


    """
    Adding methods columns
    """
    minimap2_stats = add_method_row(minimap2_stats, "minimap2")
    nanogate_positives = add_method_row(nanogate_positives, "nanogate")
    combinatorial_stats = add_method_row(combinatorial_stats, "minimap2 combinatorial")
    nanogate_time = add_method_row(nanogate_time, "nanogate")
    minimap2_constructs_table = add_method_row(minimap2_constructs_table,"minimap2")
    nanogate_constructs_table = add_method_row(nanogate_constructs_table, "nanogate")
    combinatorial_constructs_table = add_method_row(combinatorial_constructs_table,"minimap2 combinatorial")
    """
    Renaming columns
    """

    nanogate_positives = nanogate_positives.rename({"depth": "Coverage",
                                                    "constructs number": "Constructs library", "parts": "Parts number"},
                                                   errors="raise", axis="columns")

    nanogate_time = nanogate_time.rename({"Depth": "Coverage", "Constructs number": "Constructs library",
                                          "Parts": "Parts number", "Computing time": "Real Time"},
                                         errors="raise", axis="columns")

    nanogate_constructs_table = nanogate_constructs_table.rename({"depth": "Coverage",
                                                    "constructs library": "Constructs library", "parts number": "Parts number"},
                                                   errors="raise", axis="columns")

    """
    precision and  recall comparison
    """
    nanogate_positives_subtable = nanogate_positives[
        ['precision','recall', 'Constructs library', 'Coverage', "Parts number", 'Method']]
    minimap2_positives_subtable = minimap2_stats[
        ['precision','recall', 'Constructs library', 'Coverage', "Parts number", 'Method']]
    combinatorial_positives_subtable = combinatorial_stats[
        ['precision','recall','Constructs library', 'Coverage', "Parts number", 'Method']]
    true_positives_comparison = pd.concat(
        [nanogate_positives_subtable, minimap2_positives_subtable, combinatorial_positives_subtable])
    #writing output
    true_positives_comparison.to_csv(true_positives_comparison_filename)

    """
    precision and recall comparison per construct
    """
    nanogate_positives_subtable = nanogate_constructs_table[
        ['precision','recall', 'Constructs library', 'Coverage', "Parts number", 'Method']]
    minimap2_positives_subtable = minimap2_constructs_table[
        ['precision','recall', 'Constructs library', 'Coverage', "Parts number", 'Method']]
    combinatorial_positives_subtable = combinatorial_constructs_table[
        ['precision','recall','Constructs library', 'Coverage', "Parts number", 'Method']]
    true_positives_comparison = pd.concat(
        [nanogate_positives_subtable, minimap2_positives_subtable, combinatorial_positives_subtable])
    #writing output
    true_positives_comparison.to_csv(constructs_comparison_filename)


    """
    Time comparison
    """
    nanogate_time_subtable = nanogate_time[['Real Time', 'Constructs library', 'Coverage', "Parts number", 'Method']]
    minimap2_time_subtable = minimap2_stats[['Real Time', 'Constructs library', 'Coverage', "Parts number", 'Method']]
    combinatorial_time_subtable = combinatorial_stats[
        ['Real Time', 'Constructs library', 'Coverage', "Parts number", 'Method']]
    time_comparisons = pd.concat([nanogate_time_subtable, minimap2_time_subtable, combinatorial_time_subtable])
    #writing output
    time_comparisons.to_csv(time_comparison_filename)



def main():
    """
    parsing arguments
    """
    argh.dispatch_commands([comparison_table
                            ])


if __name__ == '__main__':
    main()

