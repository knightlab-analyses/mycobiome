import pandas as pd
from glob import glob
from os.path import basename

prep_files = {10632: '/qmounts/qiita_data/templates/13722_prep_10632_20210405-101218.txt',
 10639: '/qmounts/qiita_data/templates/13722_prep_10639_20210405-221712.txt',
 10648: '/qmounts/qiita_data/templates/13722_prep_10648_20210405-223449.txt',
 10649: '/qmounts/qiita_data/templates/13722_prep_10649_20210405-223544.txt',
 10650: '/qmounts/qiita_data/templates/13722_prep_10650_20210405-223628.txt',
 10651: '/qmounts/qiita_data/templates/13722_prep_10651_20210405-223710.txt',
 10652: '/qmounts/qiita_data/templates/13722_prep_10652_20210405-223753.txt',
 10653: '/qmounts/qiita_data/templates/13722_prep_10653_20210405-223831.txt',
 10654: '/qmounts/qiita_data/templates/13722_prep_10654_20210405-223915.txt',
 10655: '/qmounts/qiita_data/templates/13722_prep_10655_20210405-223959.txt',
 10656: '/qmounts/qiita_data/templates/13722_prep_10656_20210405-224101.txt',
 10820: '/qmounts/qiita_data/templates/13767_prep_10820_20210428-142955.txt',
 10821: '/qmounts/qiita_data/templates/13767_prep_10821_20210428-143047.txt',
 10822: '/qmounts/qiita_data/templates/13767_prep_10822_20210428-143129.txt',
 10823: '/qmounts/qiita_data/templates/13767_prep_10823_20210428-143445.txt',
 10824: '/qmounts/qiita_data/templates/13767_prep_10824_20210428-143524.txt',
 10825: '/qmounts/qiita_data/templates/13767_prep_10825_20210428-143604.txt',
 10826: '/qmounts/qiita_data/templates/13767_prep_10826_20210428-143652.txt',
 10827: '/qmounts/qiita_data/templates/13767_prep_10827_20210428-144015.txt',
 10828: '/qmounts/qiita_data/templates/13767_prep_10828_20210428-144102.txt',
 10829: '/qmounts/qiita_data/templates/13767_prep_10829_20210428-144400.txt',
 10830: '/qmounts/qiita_data/templates/13767_prep_10830_20210428-144440.txt',
 10831: '/qmounts/qiita_data/templates/13767_prep_10831_20210428-144535.txt',
 10833: '/qmounts/qiita_data/templates/13767_prep_10833_20210428-144942.txt',
 10834: '/qmounts/qiita_data/templates/13767_prep_10834_20210428-145027.txt',
 10835: '/qmounts/qiita_data/templates/13767_prep_10835_20210428-145208.txt',
 10836: '/qmounts/qiita_data/templates/13767_prep_10836_20210428-145504.txt',
 10837: '/qmounts/qiita_data/templates/13767_prep_10837_20210428-145545.txt',
 10838: '/qmounts/qiita_data/templates/13767_prep_10838_20210428-145805.txt',
 10839: '/qmounts/qiita_data/templates/13767_prep_10839_20210428-150151.txt',
 10840: '/qmounts/qiita_data/templates/13767_prep_10840_20210428-150226.txt',
 10841: '/qmounts/qiita_data/templates/13767_prep_10841_20210428-150258.txt',
 10842: '/qmounts/qiita_data/templates/13767_prep_10842_20210428-150600.txt',
 10843: '/qmounts/qiita_data/templates/13767_prep_10843_20210428-150648.txt',
 10844: '/qmounts/qiita_data/templates/13767_prep_10844_20210428-150837.txt',
 10845: '/qmounts/qiita_data/templates/13767_prep_10845_20210428-151005.txt',
 10846: '/qmounts/qiita_data/templates/13767_prep_10846_20210428-151038.txt',
 10847: '/qmounts/qiita_data/templates/13767_prep_10847_20210428-151115.txt',
 10848: '/qmounts/qiita_data/templates/13767_prep_10848_20210428-151356.txt',
 10849: '/qmounts/qiita_data/templates/13767_prep_10849_20210428-151432.txt',
 10850: '/qmounts/qiita_data/templates/13767_prep_10850_20210428-151506.txt',
 10851: '/qmounts/qiita_data/templates/13767_prep_10851_20210428-151721.txt',
 10852: '/qmounts/qiita_data/templates/13767_prep_10852_20210428-151751.txt',
 10853: '/qmounts/qiita_data/templates/13767_prep_10853_20210428-151824.txt',
 10854: '/qmounts/qiita_data/templates/13767_prep_10854_20210428-151858.txt',
 10855: '/qmounts/qiita_data/templates/13767_prep_10855_20210428-152032.txt',
 10856: '/qmounts/qiita_data/templates/13767_prep_10856_20210428-152118.txt',
 10857: '/qmounts/qiita_data/templates/13767_prep_10857_20210428-152153.txt',
 10858: '/qmounts/qiita_data/templates/13767_prep_10858_20210428-152225.txt',
 10859: '/qmounts/qiita_data/templates/13767_prep_10859_20210428-152300.txt',
 10860: '/qmounts/qiita_data/templates/13767_prep_10860_20210428-152332.txt',
 10861: '/qmounts/qiita_data/templates/13767_prep_10861_20210428-152357.txt',
 10862: '/qmounts/qiita_data/templates/13767_prep_10862_20210428-152422.txt',
 10863: '/qmounts/qiita_data/templates/13767_prep_10863_20210428-152450.txt',
 10889: '/qmounts/qiita_data/templates/13767_prep_10889_20210504-105940.txt'}


folders = [
    '10632-115207', '10651-115224', '10656-115229', '10824-117420', '10829-117425', 
    '10835-117432', '10840-117437', '10845-117441', '10850-117446', '10855-117451', 
    '10860-117456', '10639-115212', '10652-115227', '10820-117416', '10825-117421', 
    '10830-117426', '10836-117433', '10841-117438', '10846-117442', '10851-117447', 
    '10856-117452', '10861-117457', '10648-115221', '10653-115226', '10821-117417', 
    '10826-117422', '10831-117427', '10837-117434', '10842-117460', '10847-117443', 
    '10852-117448', '10857-117453', '10862-117458', '10649-115222', '10654-115225', 
    '10822-117418', '10827-117423', '10833-117430', '10838-117435', '10843-117439', 
    '10848-117444', '10853-117449', '10858-117454', '10863-117459', '10650-115223', 
    '10655-115228', '10823-117419', '10828-117424', '10834-117431', '10839-117436', 
    '10844-117440', '10849-117445', '10854-117450', '10859-117455', '10889-118022']
 
#
# filtered
#

taxonomy = dict()
tables = {
    'Observed_markers': [],
    'Read_counts': [],
    'Percent_observed_markers': [],
    'Total_marker_coverage': [],
    'Percent_identity': []
}
for folder in folders:
    prep_id = int(folder.split('-')[0])
    prep = pd.read_csv(prep_files[prep_id], sep='\t', dtype=str)
    prefix_to_name = prep.set_index('run_prefix')['sample_name'].to_dict()
    # making sure there are no empty files
    for htable in glob(f'Euk-{folder}/*_filtered_hits_table.txt'):
        hits = pd.read_csv(htable, sep='\t')
        if hits.shape[0] == 0:
            continue
        hits.set_index('Taxid', inplace=True)
        sname = basename(htable.replace('_filtered_hits_table.txt', '.'))
        if sname not in prefix_to_name:
            sname = basename(htable.replace('_filtered_hits_table.txt', ''))
        taxonomy.update(hits.Name.to_dict())
        for k in tables.keys():
             data = {x: float(str(y).replace('%', '')) for x, y in hits[k].to_dict().items()}
             data['sample_name'] = prefix_to_name[sname]
             tables[k].append(data)
for k, data in tables.items():
    table = pd.DataFrame(data).transpose()
    table.fillna(0, inplace=True)
    table.columns = table.loc['sample_name'].values
    table.drop('sample_name', inplace=True)
    table.index.name = '#OTU ID'
    table.to_csv(f'tables/filtered_{k}.tsv', sep='\t')

with open('tables/filtered_taxonomy.tsv', 'w') as t:
    t.write('#OTU ID\tTaxonomy\n')
    for x, y in taxonomy.items():
        t.write(f'{x}\t{y}\n')

#
# not filtered
#

taxonomy = dict()
tables = {
    'Observed_markers': [],
    'Read_counts': [],
    'Percent_observed_markers': [],
    'Total_marker_coverage': [],
    'Percent_identity': []
}
for folder in folders:
    prep_id = int(folder.split('-')[0])
    prep = pd.read_csv(prep_files[prep_id], sep='\t', dtype=str)
    prefix_to_name = prep.set_index('run_prefix')['sample_name'].to_dict()
    # making sure there are no empty files
    for htable in glob(f'Euk-{folder}/filtering/*.filtered_all_hits_table.txt'):
        hits = pd.read_csv(htable, sep='\t')
        if hits.shape[0] == 0:
            continue
        hits.set_index('Taxid', inplace=True)
        sname = basename(htable.replace('_all_hits_table.txt', '.'))
        if sname not in prefix_to_name:
            sname = basename(htable.replace('_all_hits_table.txt', ''))
        taxonomy.update(hits.Name.to_dict())
        for k in tables.keys():
            data = {x: float(str(y).replace('%', '')) for x, y in hits[k].to_dict().items()}
            data['sample_name'] = prefix_to_name[sname]
            tables[k].append(data)
for k, data in tables.items():
    table = pd.DataFrame(data).transpose()
    table.fillna(0, inplace=True)
    table.columns = table.loc['sample_name'].values
    table.drop('sample_name', inplace=True)
    table.index.name = '#OTU ID'
    table.to_csv(f'tables/not_filtered_{k}.tsv', sep='\t')

with open('tables/not_filtered_taxonomy.tsv', 'w') as t:
    t.write('#OTU ID\tTaxonomy\n')
    for x, y in taxonomy.items():
        t.write(f'{x}\t{y}\n')

