#!/usr/bin/python3.5

import re
import csv
import sys
import requests
from time import sleep
from bs4 import BeautifulSoup

v_info = sys.version_info
if v_info < (3, 5):
	raise Exception('Python version must be 3.5 at least. Your version: {}.{}'.format(v_info.major, v_info.minor))

def write_csv(filepath, headers, data):
	with open(filepath, 'w') as f:
		writer = csv.DictWriter(f, fieldnames=headers)
		writer.writeheader()
		for d in data:
			writer.writerow(d)

def parse_table(table_name, table):
	search = iter(table.select("tr"))
	data = []
	headers = [th.text.lower().replace(' ', '').replace('-', '').replace('(', '').replace(')','') for th in next(search).find_all('td')]
	for tr in search:
		data.append(dict(zip(headers, [td.text for td in tr.find_all('td')])))
	write_csv(outputs[table_name], headers, data)

if __name__ == "__main__":

	# regular expression to match with result addresses
	result_regexp = re.compile(r'https://snp-nexus\.org/test/snpnexus_[0-9]+/results\.html')

	# args
	args = sys.argv
	assert len(args) == 5
	snp_file = args[1]
	output_integrative_score_known = args[2]
	output_integrative_score_novel = args[3]
	max_retry = int(args[4]) # nb of minutes

	# file
	files = {'region_file': open(snp_file, 'rb')}

	# output file map
	outputs = {
		'novel': output_integrative_score_novel,
		'known': output_integrative_score_known
	}

	# form payload
	payload = {
	    'ensembl': 'ensembl',
		'assembly': 'hg19',
	    'email': '',
		'dataset':'',
		'query': 'batch',
		'Type': 'Chromosome',
		'Clone':'',
		'Position_clone':'',
		'a1':'',
		'a2':'',
		'strand': '1',
		'chrom': '1',
		'region_start':'',
		'region_end':'',
		'dbsnp_id':'',
		'batch_text':'',
	    'ncsnp': 'ncsnp',
	    'eigen': 'eigen',
	    'fathmm': 'fathmm',
	    'gwava': 'gwava',
	    'deepsea': 'deepsea',
	    'funseq2': 'funseq2',
	    'remm': 'remm',
	    'pipeline': 'all_var',
		'vcf': 'txt'
	}

	# For hg19 non-coding annotation, send payload and file
	url = 'https://snp-nexus.org/cgi-bin/snp/s6_nc.cgi'
	r = requests.post(url=url, data=payload, files=files)
	result = r.text

	# get link with regex
	link = None
	if result:
	    m = result_regexp.search(result)
	    if m:
	        link = m.group(0)

	# wait for results link
	html_data = None
	if link is not None:
		for i in range(max_retry):
			print("Next try: {} min...".format(i))
			sleep(60)
			r = requests.get(url=link)
			parsed_html = BeautifulSoup(r.text, "lxml")
			_table = parsed_html.body.find('table', attrs={'class':'sofT2'})
			if _table:
				html_data = BeautifulSoup(r.text, "lxml")
				break

	# # TEST PURPOSE
	# with open('./resultat_test.html') as f:
	# 	html_data = BeautifulSoup(f.read(-1), "lxml")

	# parse & csv-ishify result tables
	if html_data:
		tables = iter(html_data.find_all('table', attrs={'id':'result-table'})[-2:])
		parse_table('known', next(tables))
		parse_table('novel', next(tables))

	else:
		print('NO DATA')


	
