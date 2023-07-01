# Copyright 2022 Ashlyn Powell

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
# and associated documentation files (the "Software"), to deal in the Software without restriction, 
# including without limitation the rights to use, copy, modify, merge, publish, distribute, 
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software 
# is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial 
# portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT 
# NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import re
import sys
from sys import argv
from matplotlib import pyplot as plt, font_manager as fm
import numpy as np
import itertools

# For Plodia, uncomment lines 30-32, 66-67 (comment 64-65), 88-93

fprop = fm.FontProperties(family = 'sans serif')
hatch_pattern = "|||"

def readFile(path):
	# Takes in path to fasta file with one or two h fibroin sequences

	with open(path) as file:
		sequences = file.read().strip().split("\n")
		seqs = []
		for seq in sequences:
			if not seq.startswith(">"):
				seqs.append(seq)
	return seqs

def organizeData(sequence):
	# Split given sequence at S blocks

	sequence = sequence.replace("-","")

	# Trichoptera regex
	s_patterns = re.findall(r"S.S.S", sequence)
	i_patterns = re.split(r"S.S.S", sequence)

	# # Plodia regex
	# s_patterns = re.findall(r"SA..A", sequence)
	# i_patterns = re.split(r"SA..A", sequence)
	
	first = i_patterns.pop(0)
	a_patterns = [first]
	for i in range(len(s_patterns)):
		pattern = s_patterns[i]  + i_patterns[i]
		a_patterns.append(pattern)
	
	patterns_dict = {}
	for pattern in a_patterns:
		patterns_dict[pattern] = a_patterns.count(pattern)

	return a_patterns, patterns_dict

def lengthList(a_patterns):
	# Generate a list of lengths of the patterns

	length = []
	for item in a_patterns:
		length.append(len(item))
	return length

def checkSimilarity(keys, seq):
	# Check the number of differences for each key
	# If similar key is found, return key, else return none

	for key in keys:
		diffs = 0
		if len(key) == len(seq):
			for i in range(len(key)):
				if key[i] != seq[i]:
					diffs += 1
		if diffs <= 5 and diffs >= 1:
			if combineSeqs(key,seq).count("X") <= 5:
		# if diffs <= 10 and diffs >= 1: # Plodia
		# 	if combineSeqs(key,seq).count("X") <= 10: # Plodia 
				return [key, combineSeqs(key,seq)]
	return None

def combineSeqs(seq1, seq2):
	newSeq = []
	for i in range(len(seq1)):
		if seq1[i] != seq2[i]:
			newSeq.append("X")
		else:
			newSeq.append(seq1[i])
	return "".join(newSeq)

def colorListSingle(a_patterns):
	# Build color dictionary and list for plotting

	possible_colors = ["#482677FF", "#FDE725FF", "#453781FF", "#fdc827",
	"#440154FF", "#55C667FF", "#DCE319FF", "#fad824", "#B8DE29FF",  "#0d0887", 
	"#238A8DFF", "#350498", "#33638DFF", "#95D840FF", "#2D708EFF", 
	"#73D055FF", "#404788FF", "#90d743", "#dae319", "#1F968BFF", "#3CBB75FF"]

	# # Alternate color palette - plodia
	# possible_colors = ["#fde725","#8b0aa5","#fada24","#e97158",
	# "#0d0887","#eb5760","#febd2a","#470b6a","#9f2a63",
	# "#b83289","#5302a3","#f0f921", "#d44842","#482878",
	# "#db5c68","#a31e9a","#350498","#f48849",
	# "#6f00a8","#fba238","#cc4778","#221150", "#f3e55d"]

	colors_key = ["#17202A"]
	hatches_key = [""]

	color_dict = {a_patterns[0]:"#17202A"}
	colors = ["#17202A"]
	
	hatches_dict = {a_patterns[0]:""}
	hatches = [""]

	i = 0
	hatch = ""

	x_patterns = [a_patterns[0]]

	for item in a_patterns[1:-1]:
		if item not in color_dict:
			similars = checkSimilarity(list(color_dict.keys()), item)
			if similars != None:
				similar = similars[0]
				new = similars[1]

				val = color_dict[similar]
				color_dict.pop(similar,None)
				color_dict[new] = val

				val2 = hatches_dict[similar]
				hatches_dict.pop(similar,None)
				hatches_dict[new] = val2

				colors.append(color_dict[new])
				hatches.append(hatches_dict[new])

				x_patterns.append(new)
				x_patterns = changeSimilars(new, x_patterns)

			else:
				x_patterns.append(item)
				color_dict[item] = possible_colors[i]
				colors.append(possible_colors[i])

				hatches_dict[item] = hatch
				hatches.append(hatch)

				colors_key.append(possible_colors[i])
				if colors_key.count(possible_colors[i]) > 1:
					hatches_key.append(hatch_pattern)
				else:
					hatches_key.append("")

				i += 1
				if i == len(possible_colors):
					i = 0
					hatch = hatch_pattern
		else:
			colors.append(color_dict[item])
			hatches.append(hatches_dict[item])
			x_patterns.append(item)

		
	x_patterns.append(a_patterns[-1])
	colors_key.append("#17202A")
	hatches_key.append("")

	color_dict[a_patterns[-1]] = "#17202A"
	colors.append("#17202A")
	hatches.append("")

	return x_patterns, colors_key, hatches_key, color_dict, colors, hatches

def colorListDouble(a_patterns1, a_patterns2):
	# Build two color lists using one combined dictionary

	possible_colors = ["#482677FF", "#FDE725FF", "#453781FF", "#fdc827",
	"#440154FF", "#55C667FF", "#DCE319FF", "#fad824", "#B8DE29FF",  "#0d0887", 
	"#238A8DFF", "#350498", "#33638DFF", "#95D840FF", "#2D708EFF", 
	"#73D055FF", "#404788FF", "#90d743", "#dae319", "#1F968BFF", "#3CBB75FF"]

	colors_key = ["#17202A"]
	hatches_key = [""]

	color_dict = {a_patterns1[0]:"#17202A", a_patterns2[0]:"#17202A"}
	colors1 = ["#17202A"]
	colors2 = ["#17202A"]
	
	hatches_dict = {a_patterns1[0]:"", a_patterns2[0]:""}
	hatches1 = [""]
	hatches2 = [""]

	i = 0
	hatch = ""

	x_patterns = [a_patterns1[0]]

	for item in a_patterns1[1:-1]:
		if item not in color_dict:
			similars = checkSimilarity(list(color_dict.keys()), item)
			if similars != None:
				similar = similars[0]
				new = similars[1]

				val = color_dict[similar]
				color_dict.pop(similar,None)
				color_dict[new] = val

				val2 = hatches_dict[similar]
				hatches_dict.pop(similar,None)
				hatches_dict[new] = val2

				colors1.append(color_dict[new])
				hatches1.append(hatches_dict[new])

				x_patterns.append(new)
				x_patterns = changeSimilars(new, x_patterns)

			else:
				x_patterns.append(item)
				color_dict[item] = possible_colors[i]
				colors1.append(possible_colors[i])

				hatches_dict[item] = hatch
				hatches1.append(hatch)

				colors_key.append(possible_colors[i])
				if colors_key.count(possible_colors[i]) > 1:
					hatches_key.append(hatch_pattern)
				else:
					hatches_key.append("")

				i += 1
				if i == len(possible_colors):
					i = 0
					hatch = hatch_pattern
		else:
			colors1.append(color_dict[item])
			hatches1.append(hatches_dict[item])
			x_patterns.append(item)

	x_patterns.append(a_patterns1[-1])
	x_patterns.append(a_patterns2[0])

	for item in a_patterns2[1:-1]:
		if item not in color_dict:
			similars = checkSimilarity(list(color_dict.keys()), item)
			if similars != None:
				similar = similars[0]
				new = similars[1]

				val = color_dict[similar]
				color_dict.pop(similar,None)
				color_dict[new] = val

				val2 = hatches_dict[similar]
				hatches_dict.pop(similar,None)
				hatches_dict[new] = val2

				colors2.append(color_dict[new])
				hatches2.append(hatches_dict[new])

				x_patterns.append(new)
				x_patterns = changeSimilars(new, x_patterns)

			else:
				x_patterns.append(item)
				color_dict[item] = possible_colors[i]
				colors2.append(possible_colors[i])

				hatches_dict[item] = hatch
				hatches2.append(hatch)

				colors_key.append(possible_colors[i])
				if colors_key.count(possible_colors[i]) > 1:
					hatches_key.append(hatch_pattern)
				else:
					hatches_key.append("")

				i += 1
				if i == len(possible_colors):
					i = 0
					hatch = hatch_pattern
		else:
			colors2.append(color_dict[item])
			hatches2.append(hatches_dict[item])
			x_patterns.append(item)

		
	x_patterns.append(a_patterns2[-1])
	colors_key.append("#17202A")
	hatches_key.append("")

	color_dict[a_patterns1[-1]] = "#17202A"
	color_dict[a_patterns2[-1]] = "#17202A"
	colors1.append("#17202A")
	colors2.append("#17202A")
	hatches1.append("")
	hatches2.append("")

	return x_patterns, colors_key, hatches_key, color_dict, colors1, colors2, hatches1, hatches2

def split(in_list):
	#Split list in half

	list1 = in_list[:int((len(a_patterns)+1)/2)]
	list2 = in_list[int((len(a_patterns)+1)/2):]
	return list1, list2

def horizontalBarPlot(a_patterns, colors, lengths, hatches, path):
	# Plot one allele, horizontal
	path = "horz_" + path
	first_half, second_half = split(a_patterns)
	first_colors, second_colors = split(colors)
	first_len, second_len = split(lengths)
	first_hatch, second_hatch = split(hatches)
	
	plt.rcParams['pdf.fonttype'] = 42
	fig, ax = plt.subplots(1, 2)
	
	plot1 = ax[0]
	plot2 = ax[1]
	plot1.set_title('', fontproperties=fprop)
	plot2.set_title('', fontproperties=fprop)
	
	fig.set_figheight(8)
	y_pos1 = np.arange(len(first_half))[::-1]
	y_pos2 = np.arange(len(second_half))[::-1]

	plot1.barh(y_pos1, first_len, align='center',
		color = first_colors, hatch = first_hatch)
	plot1.set_xlim([0,max(lengths)+2])
	plot1.set_yticks([])
	plot1.set_xlabel("Number of Residues")


	plot2.barh(y_pos2, second_len, align='center',
		color = second_colors, hatch = second_hatch)
	plot2.set_xlim([0,max(lengths)+2])
	plot2.set_yticks([])
	plot2.set_xlabel("Number of Residues")


	
	# plt.show()
	plt.savefig(path)

def verticalBarPlotSingle(a_patterns, colors, lengths, path, hatches):
	# Plot one allele

	plt.rcParams['pdf.fonttype'] = 42
	fig, ax = plt.subplots()
	
	fig.set_figheight(6)
	fig.set_figwidth(35)
	x_pos = np.arange(len(a_patterns))

	ax.bar(x_pos, lengths, align='center', color = colors, hatch = hatches)
	ax.set_xticks([])
	ax.set_xlabel("Length of Allele")
	ax.set_ylabel("Number of Residues")
	ax.set_title('', fontproperties=fprop)

	# plt.show()
	plt.savefig(path)

def verticalBarPlotDouble(a_patterns1, a_patterns2, colors1, colors2, lengths1, lengths2, hatches1, hatches2, path):
	# Plot two alleles, colors coordinating

	plt.rcParams['pdf.fonttype'] = 42
	fig, ax = plt.subplots(2,1)

	plot1 = ax[0]
	plot2 = ax[1]

	plot1.set_title('', fontproperties=fprop)
	plot2.set_title('', fontproperties=fprop)
	
	fig.set_figheight(7)
	fig.set_figwidth(35)

	x_pos1 = np.arange(len(a_patterns1))
	x_pos2 = np.arange(len(a_patterns2))

	maxY = max([max(lengths1), max(lengths2)]) + 5
	maxX = max([max(x_pos1), max(x_pos2)]) + 5

	plot1.bar(x_pos1, lengths1, align='center', color = colors1, hatch = hatches1)
	plot1.set_ylim([-5,maxY])
	plot1.set_xlim([-5,maxX])
	plot1.set_xticks([])
	plot1.set_ylabel("Number of Residues")

	plot2.bar(x_pos2, lengths2, align='center', color = colors2, hatch = hatches2)
	plot2.set_ylim([-5,maxY])
	plot2.set_xlim([-5,maxX])
	plot2.set_xticks([])
	plot2.set_xlabel("Length of Allele")
	plot2.set_ylabel("Number of Residues")

	# plt.show()
	plt.savefig(path)

def outputPatterns(a_patterns, path):
	with open(path, "w") as file:
		for pattern in a_patterns:
			file.write(pattern + "\n")

def addCounts(patterns, x_patterns):
	for i, pattern in enumerate(patterns):
		pattern = f"{pattern} ({x_patterns.count(pattern)})"
		patterns[i] = pattern
	return patterns

def replace(my_list, old, new):
	for i,j in enumerate(my_list):
		if j == old:
			my_list[i] = new
	return my_list

def changeSimilars(new, x_patterns):
	for i, pattern in enumerate(x_patterns):
		similars = checkSimilarity([pattern], new)
		if similars != None:
			similar = similars[0]
			new = similars[1]
			replace(x_patterns, similar, new)
	return x_patterns

def unique(my_list):
	new_list = []
	for x in my_list:
		if x not in new_list:
			new_list.append(x)
	return new_list

def makeLegend(patterns, colors, hatches, path):

	fig, ax = plt.subplots()
	
	fig.set_figwidth(.5)
	fig.set_figheight(9)
	plt.rcParams['pdf.fonttype'] = 42

	y_pos = np.arange(len(colors))
	x_pos = list(itertools.repeat(3, len(colors)))

	ax.barh(patterns[::-1], x_pos, align='center', color=colors[::-1], hatch = hatches[::-1])

	right = ax.spines["right"]
	right.set_visible(False)
	left = ax.spines["left"]
	left.set_visible(False)
	top = ax.spines["top"]
	top.set_visible(False)
	bottom = ax.spines["bottom"]
	bottom.set_visible(False)

	ax.set_title('', fontproperties=fprop)
	ax.set_xlim([0,4.5])
	ax.set_xticks([])
	ax.tick_params(axis='y', which='major', pad=15)
	ax.yaxis.tick_right()

	plt.savefig(path, bbox_inches="tight")



if __name__ == "__main__":
	filename = argv[1]
	seqs = readFile(filename)
	if len(seqs) == 1:
		outfile = filename.split(".")[0] + "_patterns.txt"
		figure_path = filename.split(".")[0] + ".pdf"
		legend_path = filename.split(".")[0] + "_legend.pdf"
		seq = seqs[0]
		a_patterns, patterns_dict = organizeData(seq)
		x_patterns, colors_key, hatches_key, color_dict, colors, hatches = colorListSingle(a_patterns)
		patterns_unique = unique(x_patterns)
		lengths = lengthList(a_patterns)

		outputPatterns(a_patterns, outfile)
		# outputPatterns(x_patterns, "x_" + outfile)
		horizontalBarPlot(a_patterns, colors, lengths, hatches, figure_path)
		# verticalBarPlotSingle(a_patterns, colors, lengths, figure_path, hatches)
		patterns_unique = addCounts(patterns_unique, x_patterns)
		makeLegend(patterns_unique, colors_key, hatches_key, legend_path)
		
	else:
		outfile1 = filename.split(".")[0] + "_patterns1.txt"
		outfile2 = filename.split(".")[0] + "_patterns2.txt"
		figure_path = filename.split(".")[0] + ".pdf"
		legend_path = filename.split(".")[0] + "_legend.pdf"
		seq1 = seqs[0]
		seq2 = seqs[1]
		a_patterns1, patterns_dict1 = organizeData(seq1)
		a_patterns2, patterns_dict2 = organizeData(seq2)
		x_patterns, colors_key, hatches_key, color_dict, colors1, colors2, hatches1, hatches2 = colorListDouble(a_patterns1, a_patterns2, )
		lengths1 = lengthList(a_patterns1)
		lengths2 = lengthList(a_patterns2)
		patterns_unique = unique(x_patterns)
		outputPatterns(a_patterns1, outfile1)
		outputPatterns(a_patterns2, outfile2)
		verticalBarPlotDouble(a_patterns1, a_patterns2, colors1, colors2, lengths1, lengths2, hatches1, hatches2, figure_path)
	
		patterns_unique = addCounts(patterns_unique, x_patterns)
		makeLegend(patterns_unique, colors_key, hatches_key, legend_path)