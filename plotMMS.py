#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: @Hwade

import os
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

CURRENT_PATH = os.getcwd()

def get_data():
	"""
	产生数据
	"""
	readFile = open(os.path.join(CURRENT_PATH,'history.csv'),'r')
	line = readFile.readline()
	data = []
	while line:
		raw_data = line.split(',')	
		x = []
		y = []
		z = []
		for j in range(len(raw_data)/3):
			x.append(float(raw_data[j*3]))
			y.append(float(raw_data[j*3+1]))
			z.append(float(raw_data[j*3+2]))
		data.append([x,y,z])
		line = readFile.readline()
	readFile.close()
	return data

def gen_data():
	dataset = get_data()
	for data in dataset:
		yield data

def update_plot(data):
	global mmplt
	mmplt.remove()
	mmplt, = ax.plot(data[0],data[1],data[2],'r.')
	return mmplt

# Attaching 3D axis to the figure
fig = plt.figure()
ax = p3.Axes3D(fig)
mmplt, = ax.plot([0],[0],[0],'r.')
# Fifty lines of random 3-D lines
#data = get_data()[0]
# Creating fifty line objects.
# NOTE: Can't pass empty arrays into 3d version of plot()
#lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]

# Setting the axes properties
ax.set_xlim3d([0.0, 10.0])
ax.set_xlabel('X')

ax.set_ylim3d([0.0, 10.0])
ax.set_ylabel('Y')

ax.set_zlim3d([0.0, 10.0])
ax.set_zlabel('Z')

ax.set_title('Molecular Motion Simulation')

# Creating the Animation object
mms_ani = animation.FuncAnimation(fig, update_plot, gen_data, interval=10, repeat=True)

plt.show()
