#!/usr/bin/env python3
# coding=utf-8

"""
Dot_plot_like_in_BLAST. See https://github.com/shelkmike/Dot_plot_like_in_BLAST

Notes:
1) All comments, except this part of the heading comment, are in Russian. Sorry, but it's somewhat easier for me to write in Russian than in English. To understand some comment, you can use Google Translate. Names of variables are usually self-explanatory, so it is often possible to understand the meaning of a piece of code without comments. In case of trouble understanding code, ask a question at https://github.com/shelkmike/Dot_plot_like_in_BLAST/issues .
2) Throughout the code I use a Hungarian notation, which means I denote meaning of a word by using special prefixes. In particular:
s_ - variables that contain strings
n_ - variables that contain numbers
l_ - lists
d_ - dictionaries
f_ - file handlers
o_ - more complex objects
Nested data structures are denoted by several letters. For example, dl_ are dictionaries of lists and ll_ are lists of lists.

===

Эта программа нужна, чтобы делать дот-плоты похожие на те, что выдаёт онлайн-BLAST. Основные отличия:
1) dot_plot_like_in_BLAST выдаёт картинки, в которых ширина и высота пропорциональны ширине и высоте Query и Subject. В то время, как онлайн-бласт выдаёт прямоугольную картинку, в которой соотношение длин сторон всегда одно и тоже. У онлайн-бласта она прямоугольная даже когда последовательности одинаковой длины. 
2) Прямые матчи dot_plot_like_in_BLAST рисует синим, а обратные оранжевым. В то время, как онлайн-бласт все рисует одним цветом.
3) Онлайн-бласт выдаёт растровые картинки, в то время как dot_plot_like_in_BLAST выдаёт картинки в трёх форматах.
	а) Картинка в формате SVG.
	б) Картинка в формате PNG.
	в) Интерактивный файл в формате html. В нём можно приближать и удалять те или иные области, и, при необходимости, сохранить в SVG картинку только с нужной областью.

dot_plot_like_in_BLAST на вход берёт два FASTA-файла, в каждом из которых должна быть всего одна последовательность.

Для работы dot_plot_like_in_BLAST нужны:
1) BLAST+ в $PATH
2) Питоновская библиотека Plotly.

Список опций можно получить, запустив "dot_plot_like_in_BLAST.py --help".
"""

import sys
import os
import re
import shutil
import datetime
import plotly


#Сначала проверяю, все ли нужные программы доступны. Все проблемы запишу в список l_unavailable_files_and_folders, и потом напечатаю его. Если пользователь допустил ошибки ещё и в командной строке, то напечатаю оба списка проблем (недоступные программы и ошибки в командной строке) сразу.
l_unavailable_files_and_folders = []
if shutil.which("blastn") is None:
	l_unavailable_files_and_folders.append("\"blastn\" is not in $PATH")
if shutil.which("makeblastdb") is None:
	l_unavailable_files_and_folders.append("\"makeblastdb\" is not in $PATH")
			
#делаю парсинг аргументов командной строки. Можно было бы использовать argparse, но когда я делаю это без библиотек, то больше возможностей для того, чтобы сделать интерфейс таким, какой мне нравится.

s_command_line = " ".join(sys.argv) #команда, которой запущен Dot_plot_like_in_BLAST, в одну строку.
s_command_line_reduced = s_command_line #то же, что s_command_line, но после того, как я распаршу какой-нибудь аргумент, я удалю его из этой строки. Если останется какой-то нераспарсенный аргумент, значит пользователь ввёл неизвестные Dot_plot_like_in_BLAST аргументы, и нужно выдать ошибку.

#инициализирую исходные значения переменных
s_path_to_the_file_with_the_first_sequence = "" #Путь к FASTA-файлу с первой последовательностью.
s_path_to_the_file_with_the_second_sequence = "" #Путь к FASTA-файлу со второй последовательностью.
s_label_for_the_horizontal_axis = "Position in the first sequence (bp)" #Надпись, которая на графике будет размещена вдоль горизонтальной оси. Относится к первой последовательности.
s_label_for_the_vertical_axis = "Position in the second sequence (bp)" #Надпись, которая на графике будет размещена вдоль горизонтальной оси. Относится ко второй последовательности.
s_blast_tool_to_use = "blastn" #Какую программу использовать. Или "blastn", или "megablast". Думал её использовать discontiguous_megablast; но это может немного запутать пользователя потому, что в опциях самого BLAST он называется "dc_megablast"; а это название мне кажется несколько непонятным.
s_maximum_evalue = "0.001" #Максимальное разрешённое e-value. Делаю эту переменную строковой, а не числовой, потому что я не уверен, что Python будет правильно работать с числовыми переменными вроде "1e-157".
n_minimum_percent_identity = 0 #Минимальный разрешённый порог сходства, в процентах. Может быть целым или нецелым числом от 0 до 100.
n_minimum_match_length = 0 #Минимальная разрешённая длина матча.
s_additional_blast_parameters = "" #Дополнительные параметры, которые пользователь хочет передать BLAST. Например, строка вида "-word_size 7 -dust no"
n_number_of_best_matches_to_draw = 1000 #Максимальное количество матчей, которые нужно рисовать. Соответствует параметру BLAST "-max_hsps".
s_vertical_axis_direction = "bottom-up" #Ориентация диаграммы. Находится ли старт осей координат слева снизу ("bottom-up"), или диаграмма перевёрнутая и старт осей координат находится слева сверху ("top-down").
n_number_of_cpu_threads_to_use = 10 #Количество потоков для выравнивания.
s_should_diagram_be_square = "no" #Нужно ли делать так, чтобы соотношение сторон диаграммы было квадратным.
n_line_width = 1 #Толщина лини матчей, в пикселях.
n_font_size = 15 #Размер шрифта. Это размер шрифта как для подписей осей, так и для чисел, расположенных вдоль осей.
s_horizontal_tick_distance = "auto" #Расстояние между засечками по горизонтальной оси. Либо "auto" (тогда Plotly выбирает его автоматически), либо число, выраженное в количестве нуклеотидов. Если это число, то, для простоты, не конвертирую его в int, а так и оставляю строкой.
s_vertical_tick_distance = "auto" #Расстояние между засечками по вертикальной оси. Либо "auto" (тогда Plotly выбирает его автоматически), либо число, выраженное в количестве нуклеотидов. Если это число, то, для простоты, не конвертирую его в int, а так и оставляю строкой.
s_path_to_the_output_folder = "./Dot_plot_like_in_BLAST__results" #Путь к выходной папке.


s_version_of_Dot_plot_like_in_BLAST = "1.4"


l_errors_in_command_line = [] #список ошибок в командной строке. Если пользователь совершил много ошибок, то Dot_plot_like_in_BLAST напишет про них все, а не только про первую встреченную.

#если нет ни одного аргумента командной строки, или есть аргумент командной строки --help, то печатаю хелп
if (len(sys.argv) == 1) or re.search(r"\s\-\-help", s_command_line):
	print("""\nDot_plot_like_in_BLAST, a tool to make dot plots.\n
Mandatory options:
1) --file_with_the_first_sequence        Path to the FASTA file with the first sequence. This file should contain only one sequence.
2) --file_with_the_second_sequence        Path to the FASTA file with the second sequence. This file should contain only one sequence.

Alignment options:
3) --threads        Number of threads. The default value is "10".
4) --blast_tool_to_use        The BLAST tool to use for alignment. Possible values are "blastn" or "megablast". The default value is "blastn".
5) --maximum_evalue        Maximum allowed e-value. The default value is "1e-3".
6) --minimum_percent_identity        Minimum allowed percent identitity. It should be a number from 0 to 100. The default value is "0".
7) --minimum_match_length        Minimum allowed match length. The default value is "0".
8) --additional_blast_parameters        Additional parameters to pass to BLAST. Should be provided in square brackets. For example: "--additional_blast_parameters [-word_size 7 -dust no]".

Diagram options:
9) --label_for_the_horizontal_axis        Label for the horizontal axis that corresponds to the first sequence. The default value is "Position in the first sequence (bp)". Should be provided in quotes.
10) --label_for_the_vertical_axis        Label for the vertical axis that corresponds to the second sequence. The default value is "Position in the second sequence (bp)". Should be provided in quotes.
11) --vertial_axis_direction        Whether the origin of the coordinate system in the diagram should be in the left bottom (the value must be "bottom-up") or in the left top (the value must be "top-down"). The default value is "bottom-up".
12) --number_of_best_matches_to_draw        If the number of BLAST matches is larger than this value, Dot_plot_like_in_BLAST will draw only this number of matches with the smallest e-values. The default value is "1000".
13) --square_diagram        By default, sides of the diagram are proportional to sequences' lengths ("--square_diagram no"). To make the diagram square, use "--square_diagram yes". When you align sequences with very different lengths, square diagrams are easier to view.
14) --line_width        width of lines that represent matches in the diagram, in pixels. The default value is "1".
15) --font_size        font size for the diagram, in pixels. The default value is "15".
16) --horizontal_tick_distance        distance between ticks on the horizontal axis. Should be either "auto" (to choose automatically), or a number in base pairs. The default value is "auto".
17) --vertical_tick_distance        distance between ticks on the vertical axis. Should be either "auto" (to choose automatically), or a number in base pairs. The default value is "auto".

Miscellaneous options:
18) --output_folder        Path to the output folder. The default value is "./Dot_plot_like_in_BLAST__results".

Informational options:
19) --help        Print this help.
20) --version        Print the version of Dot_plot_like_in_BLAST.

Example 1 (simple):
python3 dot_plot_like_in_BLAST.py --file_with_the_first_sequence sequence_1.fasta --file_with_the_second_sequence sequence_2.fasta

Example 2 (advanced):
python3 dot_plot_like_in_BLAST.py --file_with_the_first_sequence sequence_1.fasta --file_with_the_second_sequence sequence_2.fasta --label_for_the_horizontal_axis "position in rps12, Homo sapiens (bp)" --label_for_the_vertical_axis "position in rps12, Mus musculus (bp)" --threads 20

Example 3 (more advanced. Essentially, this switches off filtering by e-value and instead filters by minimum percent identity and minimum match length):
python3 dot_plot_like_in_BLAST.py --file_with_the_first_sequence sequence_1.fasta --file_with_the_second_sequence sequence_2.fasta --label_for_the_horizontal_axis "position in rps12, Homo sapiens (bp)" --label_for_the_vertical_axis "position in rps12, Mus musculus (bp)" --threads 20 --minimum_percent_identity 80 --minimum_match_length 100 --maximum_evalue 1000000
	""")
	sys.exit()
	
#смотрю, запросил ли пользователь версию Dot_plot_like_in_BLAST
if (len(sys.argv) == 1) or re.search(r"\s\-\-version", s_command_line):
	print("Dot_plot_like_in_BLAST " + s_version_of_Dot_plot_like_in_BLAST)
	sys.exit()

#смотрю, указал ли пользователь в командной строке дополнительные параметры BLAST.
o_regular_expression_results = re.search(r" --additional_blast_parameters \[(.*?)\]", s_command_line_reduced)
if o_regular_expression_results:
	s_additional_blast_parameters = o_regular_expression_results.group(1)
	
	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)
	
	#Проверяю, указал ли пользователь какие-то параметры BLAST, которые я и так использую с BLAST. В случае, если пользователь их указал, это вызовет ошибку BLAST.
	if (re.search(r"^\-task", s_additional_blast_parameters) or re.search(r"\s\-task", s_additional_blast_parameters)):
		l_errors_in_command_line.append("You have provided the option \"-task\" through \"--additional_blast_parameters\". However, Dot_plot_like_in_BLAST provides this option to BLAST by itself")
	if (re.search(r"^\-query", s_additional_blast_parameters) or re.search(r"\s\-query", s_additional_blast_parameters)):
		l_errors_in_command_line.append("You have provided the option \"-query\" through \"--additional_blast_parameters\". However, Dot_plot_like_in_BLAST provides this option to BLAST by itself")
	if (re.search(r"^\-db", s_additional_blast_parameters) or re.search(r"\s\-db", s_additional_blast_parameters)):
		l_errors_in_command_line.append("You have provided the option \"-db\" through \"--additional_blast_parameters\". However, Dot_plot_like_in_BLAST provides this option to BLAST by itself")
	if (re.search(r"^\-out", s_additional_blast_parameters) or re.search(r"\s\-out", s_additional_blast_parameters)):
		l_errors_in_command_line.append("You have provided the option \"-out\" through \"--additional_blast_parameters\". However, Dot_plot_like_in_BLAST provides this option to BLAST by itself")
	if (re.search(r"^\-evalue", s_additional_blast_parameters) or re.search(r"\s\-evalue", s_additional_blast_parameters)):
		l_errors_in_command_line.append("You have provided the option \"-evalue\" through \"--additional_blast_parameters\". However, Dot_plot_like_in_BLAST provides this option to BLAST by itself")
	if (re.search(r"^\-line_length", s_additional_blast_parameters) or re.search(r"\s\-line_length", s_additional_blast_parameters)):
		l_errors_in_command_line.append("You have provided the option \"-line_length\" through \"--additional_blast_parameters\". However, Dot_plot_like_in_BLAST provides this option to BLAST by itself")
	if (re.search(r"^\-perc_identity", s_additional_blast_parameters) or re.search(r"\s\-perc_identity", s_additional_blast_parameters)):
		l_errors_in_command_line.append("You have provided the option \"-perc_identity\" through \"--additional_blast_parameters\". However, Dot_plot_like_in_BLAST provides this option to BLAST by itself")
	if (re.search(r"^\-qcov_hsp_perc", s_additional_blast_parameters) or re.search(r"\s\-qcov_hsp_perc", s_additional_blast_parameters)):
		l_errors_in_command_line.append("You have provided the option \"-qcov_hsp_perc\" through \"--additional_blast_parameters\". However, Dot_plot_like_in_BLAST provides this option to BLAST by itself")
	if (re.search(r"^\-max_target_seqs", s_additional_blast_parameters) or re.search(r"\s\-max_target_seqs", s_additional_blast_parameters)):
		l_errors_in_command_line.append("You have provided the option \"-max_target_seqs\" through \"--additional_blast_parameters\". However, Dot_plot_like_in_BLAST provides this option to BLAST by itself")
	if (re.search(r"^\-max_hsps", s_additional_blast_parameters) or re.search(r"\s\-max_hsps", s_additional_blast_parameters)):
		l_errors_in_command_line.append("You have provided the option \"-max_hsps\" through \"--additional_blast_parameters\". However, Dot_plot_like_in_BLAST provides this option to BLAST by itself")
	if (re.search(r"^\-num_threads", s_additional_blast_parameters) or re.search(r"\s\-num_threads", s_additional_blast_parameters)):
		l_errors_in_command_line.append("You have provided the option \"-num_threads\" through \"--additional_blast_parameters\". However, Dot_plot_like_in_BLAST provides this option to BLAST by itself")
	

#смотрю, дал ли пользователь файл с первой последовательностью. Также, проверяю, что файл с ней существует, и что в нём ровно одна последовательность.
o_regular_expression_results = re.search(r" --file_with_the_first_sequence (\S+)", s_command_line_reduced)
if o_regular_expression_results:
	s_path_to_the_file_with_the_first_sequence = o_regular_expression_results.group(1)
	
	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)
	
	if not os.path.isfile(s_path_to_the_file_with_the_first_sequence):
		l_errors_in_command_line.append("The file with the first sequence " + s_path_to_the_file_with_the_first_sequence + " does not exist.")
	#Если файл существует, то я проверяю, что там только одна строка, которая начинается с ">". Если таких строк ни одной, то считаю, что это не FASTA, и добавляю это в список проблем. Если таких строк больше одной, то считаю, что это FASTA-файл с несколькими последовательностями, и добавляю это в список проблем.
	else:
		n_number_of_strings_starting_with_the_less_than_character = 0 #Количество строк, начинающихся с ">".
		
		f_infile = open(s_path_to_the_file_with_the_first_sequence, "r")
		for s_line in f_infile:
			o_regular_expression_results = re.search(r"^>", s_line)
			
			if o_regular_expression_results:
				n_number_of_strings_starting_with_the_less_than_character += 1
		
		f_infile.close()
		
		if n_number_of_strings_starting_with_the_less_than_character == 0:
			l_errors_in_command_line.append("The file with the first sequence " + s_path_to_the_file_with_the_first_sequence + " is not in the FASTA format.")
		
		if n_number_of_strings_starting_with_the_less_than_character > 1:
			l_errors_in_command_line.append("The file with the first sequence " + s_path_to_the_file_with_the_first_sequence + " contains more than one sequence. Dot_plot_like_in_BLAST only accepts files with a single sequence.")
	
	

#смотрю, дал ли пользователь файл со второй последовательностью. Также, проверяю, что файл с ней существует, и что в нём ровно одна последовательность.
o_regular_expression_results = re.search(r" --file_with_the_second_sequence (\S+)", s_command_line_reduced)
if o_regular_expression_results:
	s_path_to_the_file_with_the_second_sequence = o_regular_expression_results.group(1)
	
	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)
	
	if not os.path.isfile(s_path_to_the_file_with_the_second_sequence):
		l_errors_in_command_line.append("The file with the second sequence " + s_path_to_the_file_with_the_second_sequence + " does not exist.")
	#Если файл существует, то я проверяю, что там только одна строка, которая начинается с ">". Если таких строк ни одной, то считаю, что это не FASTA, и добавляю это в список проблем. Если таких строк больше одной, то считаю, что это FASTA-файл с несколькими последовательностями, и добавляю это в список проблем.
	else:
		n_number_of_strings_starting_with_the_less_than_character = 0 #Количество строк, начинающихся с ">".
		
		f_infile = open(s_path_to_the_file_with_the_second_sequence, "r")
		for s_line in f_infile:
			o_regular_expression_results = re.search(r"^>", s_line)
			
			if o_regular_expression_results:
				n_number_of_strings_starting_with_the_less_than_character += 1
		
		f_infile.close()
		
		if n_number_of_strings_starting_with_the_less_than_character == 0:
			l_errors_in_command_line.append("The file with the second sequence " + s_path_to_the_file_with_the_second_sequence + " is not in the FASTA format.")
		
		if n_number_of_strings_starting_with_the_less_than_character > 1:
			l_errors_in_command_line.append("The file with the second sequence " + s_path_to_the_file_with_the_second_sequence + " contains more than one sequence. Dot_plot_like_in_BLAST only accepts files with a single sequence.")
	

#смотрю, указал ли пользователь надпись, которая должна быть вдоль горизонтальной оси. Считаю, что аргументом опции --label_for_the_horizontal_axis является всё, что идёт до начала следующей опции (то есть, до \s+\-\-) или до конца строки (то есть, до \s*$)
o_regular_expression_results = re.search(r"( --label_for_the_horizontal_axis )(.*?)(\s+\-\-|\s*$)", s_command_line_reduced)
if o_regular_expression_results:
	s_label_for_the_horizontal_axis = o_regular_expression_results.group(2)

	s_string_to_remove = re.escape(o_regular_expression_results.group(1) + o_regular_expression_results.group(2))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#смотрю, указал ли пользователь надпись, которая должна быть вдоль вертикальной оси. Считаю, что аргументом опции --label_for_the_vertical_axis является всё, что идёт до начала следующей опции (то есть, до \s+\-\-) или до конца строки (то есть, до \s*$)
o_regular_expression_results = re.search(r"( --label_for_the_vertical_axis )(.*?)(\s+\-\-|\s*$)", s_command_line_reduced)
if o_regular_expression_results:
	s_label_for_the_vertical_axis = o_regular_expression_results.group(2)

	s_string_to_remove = re.escape(o_regular_expression_results.group(1) + o_regular_expression_results.group(2))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#смотрю, указал ли пользователь название алгоритма, который нужно использовать ("blastn" или "megablast"). Также, проверяю, что пользователь указал именно одно из этих двух значений.
o_regular_expression_results = re.search(r" --blast_tool_to_use (\S+)", s_command_line_reduced)
if o_regular_expression_results:
	s_blast_tool_to_use = o_regular_expression_results.group(1)
	
	if (s_blast_tool_to_use != "blastn") and (s_blast_tool_to_use != "megablast"):
		l_errors_in_command_line.append("You indicated --blast_tool_to_use " + s_blast_tool_to_use + " but the allowed valus are \"blastn\" or \"megablast\"")

	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#смотрю, указал ли пользователь максимальное e-value.
o_regular_expression_results = re.search(r" --maximum_evalue ([\d\.eE\+\-]+)", s_command_line_reduced)
if o_regular_expression_results:
	s_maximum_evalue = o_regular_expression_results.group(1)

	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#смотрю, указал ли пользователь минимальный процент сходства.
o_regular_expression_results = re.search(r" --minimum_percent_identity ([\d\.eE\+\-]+)", s_command_line_reduced)
if o_regular_expression_results:
	n_minimum_percent_identity = float(o_regular_expression_results.group(1))

	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#смотрю, указал ли пользователь минимальную длину матча.
o_regular_expression_results = re.search(r" --minimum_match_length ([\d\.eE\+\-]+)", s_command_line_reduced)
if o_regular_expression_results:
	n_minimum_match_length = int(o_regular_expression_results.group(1))

	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#смотрю, указал ли пользователь максимальное количество матчей, которые нужно рисовать.
o_regular_expression_results = re.search(r" --number_of_best_matches_to_draw ([\d\.eE\+\-]+)", s_command_line_reduced)
if o_regular_expression_results:
	n_number_of_best_matches_to_draw = int(o_regular_expression_results.group(1))

	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#смотрю, указал ли пользователь ориентацию диаграммы.
o_regular_expression_results = re.search(r" --vertial_axis_direction (bottom\-up|top\-down)", s_command_line_reduced)
if o_regular_expression_results:
	s_vertical_axis_direction = o_regular_expression_results.group(1)

	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#смотрю, указал ли пользователь количество потоков.
o_regular_expression_results = re.search(r" --threads ([\d\.eE\+\-]+)", s_command_line_reduced)
if o_regular_expression_results:
	n_number_of_cpu_threads_to_use = int(o_regular_expression_results.group(1))

	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#смотрю, указал ли пользователь, квадратная ли должна быть диаграмма, или нет.
o_regular_expression_results = re.search(r" --square_diagram (yes|no)", s_command_line_reduced)
if o_regular_expression_results:
	s_should_diagram_be_square = o_regular_expression_results.group(1)

	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#смотрю, указал ли пользователь толщину линии на графике.
o_regular_expression_results = re.search(r" --line_width ([\d\.eE\+\-]+)", s_command_line_reduced)
if o_regular_expression_results:
	n_line_width = float(o_regular_expression_results.group(1))

	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#смотрю, указал ли пользователь размер шрифта на графике.
o_regular_expression_results = re.search(r" --font_size ([\d\.eE\+\-]+)", s_command_line_reduced)
if o_regular_expression_results:
	n_font_size = int(o_regular_expression_results.group(1))

	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#смотрю, указал ли пользователь расстояние между засечками по горизонтальной оси.
o_regular_expression_results = re.search(r" --horizontal_tick_distance (auto|[\d\.eE\+\-]+)", s_command_line_reduced)
if o_regular_expression_results:
	s_horizontal_tick_distance = o_regular_expression_results.group(1)

	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#смотрю, указал ли пользователь расстояние между засечками по вертикальной оси.
o_regular_expression_results = re.search(r" --vertical_tick_distance (auto|[\d\.eE\+\-]+)", s_command_line_reduced)
if o_regular_expression_results:
	s_vertical_tick_distance = o_regular_expression_results.group(1)

	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#смотрю, указал ли пользователь выходную папку
o_regular_expression_results = re.search(r" --output_folder (\S+)", s_command_line_reduced)
if o_regular_expression_results:
	s_path_to_the_output_folder = o_regular_expression_results.group(1)
	
	#Если в начале s_path_to_the_output_folder не стоит "." или "/", то, видимо, пользователь имеет в виду подпапку текущей папки, но не указал "./" в начале. В таком случае я добавлю "./" в начало, иначе могут быть проблемы. Не уверен насчёт Dot_plot_like_in_BLAST, но у calculate_AG это вызывало проблемы.
	if (not re.search(r"^\.", s_path_to_the_output_folder)) and (not re.search(r"^\/", s_path_to_the_output_folder)):
		s_path_to_the_output_folder = "./" + s_path_to_the_output_folder
	
	s_string_to_remove = re.escape(o_regular_expression_results.group(0))
	s_command_line_reduced = re.sub(s_string_to_remove, "", s_command_line_reduced, 1)

#Проверяю, что выходной папки не существует, либо она существует и пустая. В противном случае, говорю пользователю, что это ошибка. Не записываю эту ошибку в список l_errors_in_command_line , а сразу останавливаю работу, потому что если выходная папка уже существует, то в неё нельзя качать файлы BUSCO.
if os.path.exists(s_path_to_the_output_folder):
	if len(os.listdir(s_path_to_the_output_folder)) != 0:
		print("Dot_plot_like_in_BLAST has stopped because the output folder already exists and is not empty.")
		sys.exit()

os.makedirs(s_path_to_the_output_folder, exist_ok = True)

#проверяю, не ввёл ли пользователь какие-то несуществующие опции. Это я определяю по тому, что после того, как я распарсил все команды, в строке s_command_line_reduced осталось что-то, кроме названия исполняемого файла Dot_plot_like_in_BLAST.
s_command_line_reduced = re.sub(r"^.*?dot_plot_like_in_BLAST(\.py)?\s*", "", s_command_line_reduced)
if s_command_line_reduced != "":
	l_errors_in_command_line.append("You have provided some options which Dot_plot_like_in_BLAST doesn't know: " + s_command_line_reduced)
	
#проверяю, были ли недоступны какие-то программы, которые нужны Dot_plot_like_in_BLAST, и были ли ошибки в командной строке. Если были какие-то из этих проблем, то пишу об этом и завершаю работу Dot_plot_like_in_BLAST.
if len(l_unavailable_files_and_folders) != 0:
	#Если ошибка была всего одна.
	if len(l_unavailable_files_and_folders) == 1:
		print("There was an error with unavailable files:")
		print(l_unavailable_files_and_folders[0])
	#Если было больше одной ошибки.
	if len(l_unavailable_files_and_folders) > 1:
		print("There were errors with unavailable files:")
		n_error_number = 0 #порядковый номер ошибки. Считается от 1.
		for s_error_text in l_unavailable_files_and_folders:
			n_error_number += 1
			print(str(n_error_number) + ") " + l_unavailable_files_and_folders[n_error_number - 1])
	
	#Печатаю пустую строку, как разделитель
	print("")
	

if len(l_errors_in_command_line) != 0:
	#Если ошибка была всего одна.
	if len(l_errors_in_command_line) == 1:
		print("There was an error in the command line of Dot_plot_like_in_BLAST:")
		print(l_errors_in_command_line[0])
	#Если было больше одной ошибки.
	if len(l_errors_in_command_line) > 1:
		print("There were errors in the command line of Dot_plot_like_in_BLAST:")
		n_error_number = 0 #порядковый номер ошибки. Считается от 1.
		for s_error_text in l_errors_in_command_line:
			n_error_number += 1
			print(str(n_error_number) + ") " + l_errors_in_command_line[n_error_number - 1])
	
	#Печатаю пустую строку, как разделитель
	print("")

if (len(l_unavailable_files_and_folders) != 0) or (len(l_errors_in_command_line) != 0):
	#Если количество ошибок с недоступными файлами и папками и количество ошибок командной строки в сумме равно 1
	if (len(l_unavailable_files_and_folders) + len(l_errors_in_command_line)) == 1: 
		print("Dot_plot_like_in_BLAST has stopped. Please, fix this error and restart Dot_plot_like_in_BLAST.")
	#Если количество ошибок с недоступными файлами и папками и количество ошибок командной строки в сумме больше 1
	if (len(l_unavailable_files_and_folders) + len(l_errors_in_command_line)) > 1:
		print("Dot_plot_like_in_BLAST has stopped. Please, fix these errors and restart Dot_plot_like_in_BLAST.")
	
	sys.exit()

f_log = open(s_path_to_the_output_folder + "/log.txt", "w", buffering=1) #f_log это общий файл с логами Dot_plot_like_in_BLAST. buffering=1 означает, что буферизация идёт только на уровне строк.
o_current_time_and_date = datetime.datetime.now()
s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
f_log.write(s_current_time_and_date + "\n")
f_log.write("Started Dot_plot_like_in_BLAST\n\n")

f_log.write("You have run Dot_plot_like_in_BLAST of version " + s_version_of_Dot_plot_like_in_BLAST + " with the following command: " + s_command_line + "\n\n")


#Определю длину первой последовательности.
n_first_sequence_length = 0
s_first_sequence = ""
f_infile = open(s_path_to_the_file_with_the_first_sequence, "r")
for s_line in f_infile:
	#Если это не строка с заголовком, то считаю, что это строка с последовательностью
	if not re.search(r"^>", s_line):
		o_regular_expression_results = re.search(r"^(.+)", s_line)
		
		s_sequence_from_this_string = o_regular_expression_results.group(1)
		#удаляю всякие пробельные символы, в том числе символ переноса строки.
		s_sequence_from_this_string = re.sub(r"\s", "", s_sequence_from_this_string)

		s_first_sequence += s_sequence_from_this_string
f_infile.close()

n_first_sequence_length = len(s_first_sequence)

#Определю длину второй последовательности.
n_second_sequence_length = 0
s_second_sequence = ""
f_infile = open(s_path_to_the_file_with_the_second_sequence, "r")
for s_line in f_infile:
	#Если это не строка с заголовком, то считаю, что это строка с последовательностью
	if not re.search(r"^>", s_line):
		o_regular_expression_results = re.search(r"^(.+)", s_line)
		
		s_sequence_from_this_string = o_regular_expression_results.group(1)
		#удаляю всякие пробельные символы, в том числе символ переноса строки.
		s_sequence_from_this_string = re.sub(r"\s", "", s_sequence_from_this_string)

		s_second_sequence += s_sequence_from_this_string
f_infile.close()

n_second_sequence_length = len(s_second_sequence)

#Делаю базу BLAST.
os.system("makeblastdb -dbtype nucl -in " + s_path_to_the_file_with_the_second_sequence + " -out " + s_path_to_the_output_folder + "/blast_database_made_from_the_second_sequence") 

o_current_time_and_date = datetime.datetime.now()
s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
f_log.write(s_current_time_and_date + "\n")
f_log.write("Started BLAST alignment\n\n")

"""
Делаю выравнивание BLAST. 
"-max_hsps" ставлю 1000000, а не n_number_of_best_matches_to_draw. Это потому, что я ещё буду отдельно делать фильтрацию по n_minimum_match_length (сам бласт по длине матча не умеет делать фильтрацию), поэтому если я установлю -max_hsps равным n_number_of_best_matches_to_draw , то в итоге матчей на рисунке может быть меньше, чем хотел пользователь.

-qcov_hsp_perc (100 * n_minimum_match_length / n_first_sequence_length) я использую затем, чтобы ограничить размер файла, выдаваемого бластом. Всё равно длина выровнявшегося фрагмента в query всегда меньше либо равна длине выравнивания. Поэтому когда я ставлю ограничение "-qcov_hsp_perc (100 * n_minimum_match_length / n_first_sequence_length)", я не выкидываю никакие матчи, которые нужны пользователю.
"""

os.system("blastn -task " + s_blast_tool_to_use + " -query " + s_path_to_the_file_with_the_first_sequence + " -db " + s_path_to_the_output_folder + "/blast_database_made_from_the_second_sequence -out " + s_path_to_the_output_folder + "/blast_results.txt -evalue " + s_maximum_evalue + " -line_length 1000000000 -perc_identity " + str(n_minimum_percent_identity) + " -qcov_hsp_perc " + str(100 * n_minimum_match_length / n_first_sequence_length) + " -max_target_seqs 1 -max_hsps 1000000 " + s_additional_blast_parameters + " -num_threads " + str(n_number_of_cpu_threads_to_use))

o_figure_with_dotplot = plotly.graph_objects.Figure() #рисунок с дотплотом

"""
Определяю ширину и высоту рисунка. 
Если пользователь указал "--square_diagram yes", то ширина и высота будут по 550 пикселей.
Если "--square_diagram no", то минимальное из двух измерений будет 550, а максимальное будет равно 550*соотношение_длин_последовательностей.
Дополнительно, по 100 пикселей с каждой из четырёх сторон будет margins (там подписи и всякое такое).

Именно 550 пикселей использую, потому что тогда высота с рамками будет 750 пикселей. Тогда при просмотре HTML картинка целиком влезает в окно браузера на компьютере, на котором вертикальное разрешение монитора 1080 пикселей (ещё некоторое расстояние нужно для панели задач, вкладок браузера, адресной строки браузера и некоторых других элементов браузера).
"""
if s_should_diagram_be_square == "yes":
	n_diagram_width = 550
	n_diagram_height = 550
else:
	if n_first_sequence_length < n_second_sequence_length:
		n_diagram_width = 550
		n_diagram_height = int(550 * n_second_sequence_length / n_first_sequence_length)
	elif n_first_sequence_length > n_second_sequence_length:
		n_diagram_width = int(550 * n_first_sequence_length / n_second_sequence_length)
		n_diagram_height = 550
	else:
		n_diagram_width = 550
		n_diagram_height = 550

o_figure_with_dotplot.update_layout(width = n_diagram_width + 100 + 100, height = n_diagram_height + 100 + 100, margin = dict(l = 100, r = 100, t = 100, b = 100))

#Определяю координаты осей и подписи осей. Также, если нужно, инвертирую ось y. Когда ось y инвертирована, то и подписи по оси x должны быть сверху.
if s_vertical_axis_direction == "top-down":
	o_figure_with_dotplot.update_layout(xaxis = {"range" : [0, n_first_sequence_length], "side" : "top"}, xaxis_title = s_label_for_the_horizontal_axis)
	o_figure_with_dotplot.update_layout(yaxis = {"range" : [n_second_sequence_length, 0]}, yaxis_title = s_label_for_the_vertical_axis)
else:
	o_figure_with_dotplot.update_layout(xaxis = {"range" : [0, n_first_sequence_length]}, xaxis_title = s_label_for_the_horizontal_axis)
	o_figure_with_dotplot.update_layout(yaxis = {"range" : [0, n_second_sequence_length]}, yaxis_title = s_label_for_the_vertical_axis)

"""
По моему опыту, если добавлять линии по одной с помощью o_figure_with_dotplot.add_shape, то Plotly работает запредельно медленно. Поэтому я использую решение, описанное на https://stackoverflow.com/questions/70276242/adding-500-circles-in-a-plotly-graph-using-add-shape-function-takes-45-seconds . Я делаю список с линиями, а затем использую o_figure_with_dotplot.update_layout , чтобы добавить их. Линии могут идти только под углом в 45 градусов или углом в 135 градусов.
"""

o_current_time_and_date = datetime.datetime.now()
s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
f_log.write(s_current_time_and_date + "\n")
f_log.write("Finished BLAST alignment. Now parsing BLAST results\n\n")

l_lines_for_input_to_Plotly = []

n_number_of_matches_that_pass_thresholds__that_I_have_seen = 0 #Количество матчей, которые я видел, и которые при этом проходят порог по минимальной длине матча. Считается от 1. Когда это значение станет равно n_number_of_best_matches_to_draw, я выхожу из цикла.

#Иду по результатам BLAST
f_infile = open(s_path_to_the_output_folder + "/blast_results.txt", "r")
l_infile = f_infile.readlines()

#Также, просто для контроля, сделаю tsv-файл, в котором записаны координаты всех линий, которые должны быть на графике. Также, там будет столбец "strand", в котором будет значение "+" или "-".
f_coordinates_of_lines = open(s_path_to_the_output_folder + "/coordinates_of_lines.tsv", "w")
f_coordinates_of_lines.write("left coordinate in sequence 1\tright coordinate in sequence 1\tleft coordinate in sequence 2\tright coordinate in sequence 2\tstrand\n")

n_line_number = 0

while n_line_number < len(l_infile):
	"""
	Когда я делал выравнивание BLAST, я поставил -line_length 1000000000, так что строки могут быть очень длинными.
	
	Query  16631  TGTTTTATTTAAAATATTATAATATCAAAATAAATTAAA  16669
					| ||||||||||||  |||| ||  | ||||||| ||
	Sbjct  1932   AATATTATTTAAAATA--ATAA-ATTTATATAAATTTAA  1897
	"""
	
	#Если это строка с выравниванием фрагмента Query
	o_regular_expression_results = re.search(r"Query\s*(\d+)\s*([A-Za-z\-\*]+)\s*(\d+)", l_infile[n_line_number])
	if o_regular_expression_results:
		n_start_of_the_query_line = int(o_regular_expression_results.group(1))
		s_query_alignment_line = o_regular_expression_results.group(2)
		n_end_of_the_query_line = int(o_regular_expression_results.group(3))
		
		#Если длина выравнивания больше, чем n_minimum_match_length , значит этот матч проходит порог по длине. В противном случае пропускаю его.
		if len(s_query_alignment_line) >= n_minimum_match_length:
			
			n_number_of_matches_that_pass_thresholds__that_I_have_seen += 1

			#Если количество виденных мной матчей, проходяхищ порог по минимальной длине матча, уже превысило n_number_of_best_matches_to_draw элементов, то выхожу из цикла.
			if (n_number_of_matches_that_pass_thresholds__that_I_have_seen > n_number_of_best_matches_to_draw):
				break
			
			#Строка с выравниванием фрагмента Subject идёт через две строки после выравнивания фрагмента Query.
			o_regular_expression_results = re.search(r"Sbjct\s*(\d+)\s*([A-Za-z\-\*]+)\s*(\d+)", l_infile[n_line_number + 2])
			if o_regular_expression_results:
				n_start_of_the_subject_line = int(o_regular_expression_results.group(1))
				s_subject_alignment_line = o_regular_expression_results.group(2)
				n_end_of_the_subject_line = int(o_regular_expression_results.group(3))
				
				s_match_type = "" #тип матча. "direct", если прямой; "reverse complement", если обратно-комплементарный.
				
				if n_end_of_the_subject_line >= n_start_of_the_subject_line:
					s_match_type = "direct"
				else:
					s_match_type = "reverse complement"
				
				s_am_I_currently_drawing_a_line = "no" #Нахожусь ли я сейчас в процессе рисования линии. То есть, иду ли внутри области, в которой у обоих последовательностей не гэпы. Если я сейчас в процессе рисования линии, то значение этой переменной "yes".
				
				#Теперь иду по последовательности матча, и добавляю в l_lines_for_input_to_Plotly линии. n_position_in_the_alignment_line это позиция в строке выравнивания; она считается от 1.
				
				n_number_of_seen_gaps_in_the_first_sequence = 0 #сколько гэпов в первой последовательности я уже видел
				n_number_of_seen_gaps_in_the_second_sequence = 0 #сколько гэпов во второй последовательности я уже видел
				
				for n_position_in_the_alignment_line in range(1, len(s_query_alignment_line) + 1):
					#Смотрю, правда ли в обеих последовательностях в этой колонке нет гэпа.
					if (s_query_alignment_line[n_position_in_the_alignment_line - 1] != "-") and (s_subject_alignment_line[n_position_in_the_alignment_line - 1] != "-"):
						
						#Если я был не в процессе рисования линии, то начинаю рисовать линию.
						if s_am_I_currently_drawing_a_line == "no":
							if s_match_type == "direct":
								n_match_left_coordinate_in_the_first_sequence = n_start_of_the_query_line + n_position_in_the_alignment_line - 1 - n_number_of_seen_gaps_in_the_first_sequence
								n_match_left_coordinate_in_the_second_sequence = n_start_of_the_subject_line + n_position_in_the_alignment_line - 1 - n_number_of_seen_gaps_in_the_second_sequence
							if s_match_type == "reverse complement":
								n_match_left_coordinate_in_the_first_sequence = n_start_of_the_query_line + n_position_in_the_alignment_line - 1 - n_number_of_seen_gaps_in_the_first_sequence
								n_match_right_coordinate_in_the_second_sequence = n_start_of_the_subject_line - n_position_in_the_alignment_line + 1 + n_number_of_seen_gaps_in_the_second_sequence
							
							s_am_I_currently_drawing_a_line = "yes"
					
					#Если я дошёл до гэпа
					else:
						#Если гэп был в первой последовательности
						if s_query_alignment_line[n_position_in_the_alignment_line - 1] == "-":
							n_number_of_seen_gaps_in_the_first_sequence += 1
						
						#Если гэп был во второй последовательности
						if s_subject_alignment_line[n_position_in_the_alignment_line - 1] == "-":
							n_number_of_seen_gaps_in_the_second_sequence += 1
					
						#Если я был не в процессе рисования линии, то ничего делать не нужно.
						if s_am_I_currently_drawing_a_line == "no":	
							pass
						#Если я был в процессе рисования линии, то нужно завершить процесс рисования линии.
						else:
							if s_match_type == "direct":
								n_match_right_coordinate_in_the_first_sequence = n_start_of_the_query_line + n_position_in_the_alignment_line - 1 - n_number_of_seen_gaps_in_the_first_sequence
								n_match_right_coordinate_in_the_second_sequence = n_start_of_the_subject_line + n_position_in_the_alignment_line - 1 - n_number_of_seen_gaps_in_the_second_sequence
								
								#Вычитаю единицу из координаты, чтобы не учитывать сам символ гэпа. В той последовательности, в которой есть гэп, я единицу и так вычел (потому что выше я для неё увеличил значение переменной n_number_of_seen_gaps_in_... и эту переменную я вычитал), так что мне нужно вычесть единицу только из координаты другой последовательности.
								if s_query_alignment_line[n_position_in_the_alignment_line - 1] == "-": #Если гэп в первой последовательности
									n_match_right_coordinate_in_the_second_sequence -= 1
								if s_subject_alignment_line[n_position_in_the_alignment_line - 1] == "-": #Если гэп во второй последовательности
									n_match_right_coordinate_in_the_first_sequence -= 1
								
								l_lines_for_input_to_Plotly.append(dict(type="line", x0 = n_match_left_coordinate_in_the_first_sequence, y0 = n_match_left_coordinate_in_the_second_sequence, x1 = n_match_right_coordinate_in_the_first_sequence, y1 = n_match_right_coordinate_in_the_second_sequence, line = dict(color = "#4472C4", width = n_line_width)))
								
								f_coordinates_of_lines.write(str(n_match_left_coordinate_in_the_first_sequence) + "\t" + str(n_match_right_coordinate_in_the_first_sequence) + "\t" + str(n_match_left_coordinate_in_the_second_sequence) + "\t" + str(n_match_right_coordinate_in_the_second_sequence) + "\t+\n")
								
							if s_match_type == "reverse complement":
								n_match_right_coordinate_in_the_first_sequence = n_start_of_the_query_line + n_position_in_the_alignment_line - 1 - n_number_of_seen_gaps_in_the_first_sequence
								n_match_left_coordinate_in_the_second_sequence = n_start_of_the_subject_line - n_position_in_the_alignment_line + 1 + n_number_of_seen_gaps_in_the_second_sequence
								
								#Вычитаю единицу из координаты первой последовательности, если гэп во второй. Если гэп во второй последовательности, то его я и так учёт, потому что сверху добавляю n_number_of_seen_gaps_in_the_second_sequence, которое ещё выше было увеличено на единицу.
								#Если гэп в первой последовательности, то, поскольку матч обратно-комплементарный, наоборот, прибавляю единицу к координате во второй. 
								if s_query_alignment_line[n_position_in_the_alignment_line - 1] == "-": #Если гэп в первой последовательности
									n_match_left_coordinate_in_the_second_sequence += 1
								if s_subject_alignment_line[n_position_in_the_alignment_line - 1] == "-": #Если гэп во второй последовательности
									n_match_right_coordinate_in_the_first_sequence -= 1
								
								l_lines_for_input_to_Plotly.append(dict(type="line", x0 = n_match_left_coordinate_in_the_first_sequence, y0 = n_match_right_coordinate_in_the_second_sequence, x1 = n_match_right_coordinate_in_the_first_sequence, y1 = n_match_left_coordinate_in_the_second_sequence, line = dict(color = "#ED7D31", width = n_line_width)))
								
								f_coordinates_of_lines.write(str(n_match_left_coordinate_in_the_first_sequence) + "\t" + str(n_match_right_coordinate_in_the_first_sequence) + "\t" + str(n_match_left_coordinate_in_the_second_sequence) + "\t" + str(n_match_right_coordinate_in_the_second_sequence) + "\t-\n")
								
							s_am_I_currently_drawing_a_line = "no"
				
				#Если, когда я подошёл к концу последовательности, я был в процессе рисования линии, то заканчиваю линию.
				if s_am_I_currently_drawing_a_line == "yes":	
					if s_match_type == "direct":
						n_match_right_coordinate_in_the_first_sequence = n_start_of_the_query_line + n_position_in_the_alignment_line - 1 - n_number_of_seen_gaps_in_the_first_sequence
						n_match_right_coordinate_in_the_second_sequence = n_start_of_the_subject_line + n_position_in_the_alignment_line - 1 - n_number_of_seen_gaps_in_the_second_sequence

						l_lines_for_input_to_Plotly.append(dict(type="line", x0 = n_match_left_coordinate_in_the_first_sequence, y0 = n_match_left_coordinate_in_the_second_sequence, x1 = n_match_right_coordinate_in_the_first_sequence, y1 = n_match_right_coordinate_in_the_second_sequence, line = dict(color = "#4472C4", width = n_line_width)))
						
						f_coordinates_of_lines.write(str(n_match_left_coordinate_in_the_first_sequence) + "\t" + str(n_match_right_coordinate_in_the_first_sequence) + "\t" + str(n_match_left_coordinate_in_the_second_sequence) + "\t" + str(n_match_right_coordinate_in_the_second_sequence) + "\t+\n")
						
					if s_match_type == "reverse complement":
						n_match_right_coordinate_in_the_first_sequence = n_start_of_the_query_line + n_position_in_the_alignment_line - 1 - n_number_of_seen_gaps_in_the_first_sequence
						n_match_left_coordinate_in_the_second_sequence = n_start_of_the_subject_line - n_position_in_the_alignment_line + 1 + n_number_of_seen_gaps_in_the_second_sequence

						l_lines_for_input_to_Plotly.append(dict(type="line", x0 = n_match_left_coordinate_in_the_first_sequence, y0 = n_match_right_coordinate_in_the_second_sequence, x1 = n_match_right_coordinate_in_the_first_sequence, y1 = n_match_left_coordinate_in_the_second_sequence, line = dict(color = "#ED7D31", width = n_line_width)))
						
						f_coordinates_of_lines.write(str(n_match_left_coordinate_in_the_first_sequence) + "\t" + str(n_match_right_coordinate_in_the_first_sequence) + "\t" + str(n_match_left_coordinate_in_the_second_sequence) + "\t" + str(n_match_right_coordinate_in_the_second_sequence) + "\t-\n")
						
					s_am_I_currently_drawing_a_line = "no"
			
	
	n_line_number += 1
	
	
f_infile.close()
f_coordinates_of_lines.close()

o_current_time_and_date = datetime.datetime.now()
s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
f_log.write(s_current_time_and_date + "\n")
f_log.write("Finished parsing BLAST results. Now started drawing the diagram.\n\n")

o_figure_with_dotplot.update_layout(shapes = l_lines_for_input_to_Plotly)

#Устанавливаю размер шрифта.
o_figure_with_dotplot.update_layout(font = {"size" : n_font_size})

#Устанавливаю расстояние между засечками по горизонтальной и вертикальной осям в случае, если значения s_horizontal_tick_distance и s_vertical_tick_distance оба не равны "auto".
if ((s_horizontal_tick_distance == "auto") and (s_vertical_tick_distance == "auto")):
	pass
elif (s_horizontal_tick_distance == "auto"): #Если пользователь указал, что по горизонтальной оси расстояние между засечками определяется автоматически, а по вертикальной — нет.
	o_figure_with_dotplot.update_layout(yaxis = {"dtick" : int(s_vertical_tick_distance)})
elif (s_vertical_tick_distance == "auto"): #Если пользователь указал, что по вертикальной оси расстояние между засечками определяется автоматически, а по горизонтальной — нет.
	o_figure_with_dotplot.update_layout(xaxis = {"dtick" : int(s_horizontal_tick_distance)})
else:
	o_figure_with_dotplot.update_layout(xaxis = {"dtick" : int(s_horizontal_tick_distance)}, yaxis = {"dtick" : int(s_vertical_tick_distance)})

#Делаю рисунок в формате PNG
o_figure_with_dotplot.write_image(s_path_to_the_output_folder + "/image.png")

#Делаю рисунок в формате SVG
o_figure_with_dotplot.write_image(s_path_to_the_output_folder + "/image.svg")

#Делаю рисунок в формате HTML
#"toImageButtonOptions ..." нужно, чтобы из html-файла картинки сохранялись не в png, а в svg.
o_figure_with_dotplot.write_html(s_path_to_the_output_folder + "/image.html", include_plotlyjs = "cdn", config = {'toImageButtonOptions': {'format': 'svg', 'filename': 'image', 'height': n_diagram_height, 'width': n_diagram_height, 'scale': 1}})


o_current_time_and_date = datetime.datetime.now()
s_current_time_and_date = o_current_time_and_date.strftime("%H:%M:%S %Y-%m-%d")
f_log.write(s_current_time_and_date + "\n")
f_log.write("Finished drawing the diagram.\n\nDot_plot_like_in_BLAST has finished its work.")
