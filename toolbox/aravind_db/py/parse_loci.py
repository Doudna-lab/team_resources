# == Native Modules
import time
from http.client import RemoteDisconnected
import requests
import re
# == Installed Modules
import pandas as pd
from bs4 import BeautifulSoup as bs
from urllib.parse import urljoin
from Bio import Entrez
import urllib.request
import urllib.error
import urllib3
# == Project Modules


def recursive_scrape(source_ftp_link):
	link_list = []
	subdir_list = []
	try:
		source_ftp_link.endswith('/')
	except TypeError:
		return link_list, subdir_list

	if source_ftp_link.endswith('/'):
		source_ftp_link = source_ftp_link[:-1]
		page = requests.get(source_ftp_link)
		if page.status_code != 200:
			raise "Page did not return 200 status code"
		# Create the soup object for easy navigation to various tags.
		soup = bs(page.content, 'html.parser')
		# links = soup.find_all('a', href=True)  # Reading all the table tags
		for link in soup.find_all('a', href=True):  # Reading all the table tags:

			if link.find_all(string="Parent Directory"):
				continue
			# Construct the absolute URL if the href link is relative
			absolute_url = source_ftp_link + "/" + link['href']
			if absolute_url.endswith("/"):
				subdir_list.append(absolute_url)
			elif re.search(r'\.html?$', absolute_url):
				if re.search("index.html", absolute_url):
					continue
				link_list.append(absolute_url)
	return link_list, subdir_list


def surface_scrape(source_ftp_links, root_url, scraping_iterations):
	html_list = []
	subfolders = []
	for link in source_ftp_links:
		absolute_url = root_url + "/" + link['href']
		if re.search(link['href'], root_url):
			continue
		loop_html_list, loop_subfolders = recursive_scrape(absolute_url)
		html_list.extend(loop_html_list)
		subfolders.extend(loop_subfolders)

	for _ in range(0, scraping_iterations):
		for folder in subfolders:
			loop_html_list, loop_subfolders = recursive_scrape(folder)
			html_list.extend(loop_html_list)
			subfolders.extend(loop_subfolders)

	return html_list, subfolders


def html_pre_to_df(soup_object):
	pre = soup_object.select_one('pre').text
	results = []
	for line in pre.split('\n')[1:-1]:
		if '--' not in line:
			row = re.split(r'\s+', line.strip())
			if len(row) > 1:
				print(row)
				results.append(row)
	print(f"Sending raw list of lines to split. Total length: {len(results)}")
	split_tables = split_html_table(results)
	return split_tables


def split_html_table(rows_list):
	tables_dictionary = {}
	table_closed = False
	previous_lower_limit_row_idx = int(0)
	upper_limit_table_idx = int(0)
	table_header_row = int(0)
	upper_id_anchor_row = int(-1)
	max_table_interval_window = int(3)
	for row_idx in range(0, len(rows_list)):
		row = "\t".join(rows_list[row_idx])
		# Try to anchor entry-based table on the presence of GI at the start
		# BUT THERE ARE ENTRIES ASSIGNED AS GI THAT START WITH LETTERS
		# TODO: ADDRESS THIS
		if re.search(r"^\d{5,}", row):
			# == Everything in the middle of any table
			if upper_id_anchor_row > 0:
				# == Always assume this could be the last entry of the table
				previous_lower_limit_row_idx = row_idx

			if table_closed:
				table_closed = False
				table_header_row = max(row_idx - 1, 0)

			# == Process information of the 1st entry in the 1st table
			if upper_limit_table_idx == 0:
				upper_id_anchor_row = row_idx
				table_header_row = max(row_idx - 1, 0)
				upper_limit_table_idx = max(row_idx - max_table_interval_window, 0)
		elif not re.search(r"^\d{5,}", row):
			# == Process information from 2nd table onwards
			if previous_lower_limit_row_idx > 0:
				table_closed = True
				upper_limit_table_idx = row_idx
				tentative_table_name_rows = rows_list[table_header_row - 1, upper_limit_table_idx]
				table_title = guess_table_name(tentative_table_name_rows)
				# == Capture table information to store them
				table_title = "_".join(rows_list[table_header_row - 1])
				table_content = list(rows_list[table_header_row:previous_lower_limit_row_idx])
				current_table_df = pd.DataFrame(table_content)
				current_table_df.columns = current_table_df.iloc[0, :]
				current_table_df = current_table_df.drop(0)
				# print(f"Table title: {table_title} of length {len(current_table_df)}")
				tables_dictionary.setdefault(table_title, current_table_df)

	return tables_dictionary


def guess_table_name(rows_list):

	pass


def cross_db_search(query_list, db_from, dbto):
	progress = 0
	source2target = {}
	not_found_list = []
	full_record = ''
	for hit in query_list:
		progress += 1
		dup_check = []
		uid = ''
		uid_list = []
		max_retries = 50
		search_record = {}
		# Standardize protein identifiers to NCBI UIDs through ESearch
		# Introduce a delay of 1 second before making the request

		# Attempt the search with retries
		for hit in query_list:
			try:
				handle = Entrez.esearch(db=f"{db_from}", term=f"{hit}",  idtype="text", retmax=5)
				search_record = Entrez.read(handle)
				print("Processed Esearch request")
				try:
					uid_list = search_record['IdList']

				except IndexError:
					print("No records found for")
					continue
				handle.close()
			except urllib.error.HTTPError as e:
				if e.code == 429:  # HTTP 429: Too Many Requests
					print(f"Received HTTP 429 error. Retrying in 10 seconds...")
					time.sleep(10)
				else:
					print(f"Error")
					continue  # Re-raise other HTTP errors
			except RuntimeError:
				print('Something went wrong. Retrying in')
				continue
			print("process uid list of length ", len(uid_list))
			# Loop through databases (found in config) and grab Nuccore UIDs
			if len(uid_list) >= 1:
				print("INSIDE LOOP")
				for uid in uid_list:
					if uid in set(dup_check):
						continue
					print(f"Process Elink routine: {uid}")
					for attempt in range(1, max_retries):
						try:
							link_list, loop_nuc_acc, not_found_hit, full_record = elink_routine(db_from, dbto, uid)
						except RemoteDisconnected:
							print(f"HTTP Issue - Retrying n {attempt}...")
							if attempt == max_retries:
								print("Reached Maximum number of retries")
								break
							continue
					if not_found_hit:
						not_found_list.append(not_found_hit)
						continue
					if link_list:
						dup_check.append(uid)
						source2target.setdefault(loop_nuc_acc, link_list)

	return source2target, list(set(not_found_list)), full_record


def elink_routine(dbfrom, dbto, hit_uid):
	dup_check = []
	not_found = ""
	linked = ""
	link_record = ""
	handle = None
	server_attempts = 0
	try:
		handle = Entrez.elink(dbfrom=f"{dbfrom}", db=dbto, id=f"{hit_uid}", idtype="uid")
	except urllib.error.HTTPError as err:
		if err.code == 500:
			print(f'An internal server error occurred while handling the accession {hit_uid}')
			not_found = hit_uid
			return linked, hit_uid, not_found
	try:
		link_record = Entrez.read(handle)
	except RuntimeError:
		not_found = hit_uid
	if link_record:
		try:
			linked = []
			# linked = link_record[0]['LinkSetDb'][0]['Link'][0]['Id']
			for link_idx in range(len(link_record[0]['LinkSetDb'][0]['Link'])):
				link_id = link_record[0]['LinkSetDb'][0]['Link'][link_idx]['Id']
				linked.append(link_id)
				if link_id not in dup_check:
					dup_check.append(link_id)
		except (IndexError, KeyError):
			not_found = hit_uid
	handle.close()
	return linked, hit_uid, not_found, link_record


def main():
	entrez_login = 'daniel.bellieny@gmail.com'
	Entrez.email = entrez_login
	# url = 'https://ftp.ncbi.nlm.nih.gov/pub/aravind'
	url = 'https://ftp.ncbi.nlm.nih.gov/pub/aravind/UB/prok_ub_supplement_file_1.html'

	page = requests.get(url)
	if page.status_code != 200:
		raise "Page did not return 200 status code"

	# === Parse FTP and underlying HTML files
	# Create the soup object for easy navigation to various tags.
	soup = bs(page.content, 'html.parser')
	# Upon inspection of the page — data to be extracted is in id='production-data-table'
	links = soup.find_all('a', href=True)  # Reading all the table tags
	# Iterate through FTP page multiple times and retrieve HTML files
	html_list, subfolder_list = surface_scrape(links, url, 5)

	# === Import HTMLs to Tables
	table_data = pd.DataFrame()

	dataframes = []
	for html_string in html_list:
		html = requests.get(html_string)
		soup = bs(html.content, 'lxml')
		tables = soup.select('table')  # Reading all the table tags

		if tables:
			page_dataframes = []  # List to store Pandas DataFrames for each table on this page
			for table in tables:
				table_head = table.find('thead')  # Get the head tag
				header_rows = table_head.find_all('tr')  # Get the tr
				# Convert each table to a Pandas DataFrame
				df = pd.read_html(str(table))
				page_dataframes.extend(df)  # Append each DataFrame to the list
			dataframes.append(page_dataframes)  # Append the list of DataFrames for this page

			print(f"Tables Found {len(tables)}")
			print(html_string)
			break

	a, b, c = cross_db_search(['Aravind L.'], 'pubmed', 'nuccore')
	pass


if __name__ == "__main__":
	main()



