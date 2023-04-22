import pandas as pd
import argh









def summary_table(output_summary_table: 'output filename',
                  *run_tables: 'run tables',
                  ):
    """
    Method to cooncat run tables
    """
    tables = []
    for filename in run_tables:
        table = pd.read_csv(filename)
        tables.append(table)
    joined_table = pd.concat(tables)
    joined_table.to_csv(output_summary_table)





def main():
    """
    parsing arguments
    """
    argh.dispatch_commands([summary_table
                            ])


if __name__ == '__main__':
    main()
