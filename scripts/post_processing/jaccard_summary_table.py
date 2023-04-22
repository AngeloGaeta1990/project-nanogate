import pandas as pd
import argh









def summary_jaccard_table(output_jaccard_table: 'output filename',
                  *jaccard_tables: 'jaccard tables',
                  ):
    """
    Method to cooncat run tables
    """
    tables = []
    for filename in jaccard_tables:
        table = pd.read_csv(filename)
        tables.append(table)
    joined_table = pd.concat(tables)
    joined_table.to_csv(output_jaccard_table)





def main():
    """
    parsing arguments
    """
    argh.dispatch_commands([summary_jaccard_table
                            ])

if __name__ == '__main__':
    main()
