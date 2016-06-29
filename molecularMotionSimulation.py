#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: @Hwade

import random 
import numpy
import math
import time
import csv
import os
from multiprocessing import Process, Lock, Queue

MOLECULAR_NUM = 1000
DEMENTION = 3
SPACE_LEN = 10
M = 1
PI = math.pi
dt = 0.0001

LOCK = Lock()
MUTEX = Lock()
QUEUE = Queue()
DATA_QUEUE = Queue()

lockCount = 2
process_pool = []

class MolecularMotionSimulator():
	''' 分子运动仿真 '''		

	def __init__(self):
		self.pos = [0.0]*MOLECULAR_NUM*DEMENTION
		self.vel = [0.0]*MOLECULAR_NUM*DEMENTION
		self.acc = [0.0]*MOLECULAR_NUM*DEMENTION
		# 保存新的加速度
		self.acc_dt = [0.0]*MOLECULAR_NUM*DEMENTION

		for i in range(DEMENTION):
			for j in range(MOLECULAR_NUM):
				self.pos[i+j*DEMENTION] = random.uniform(0,10)

	def direction(self, vetor):
		# 返回三个方向的余弦值
		sqt = numpy.sqrt(vetor[0]**2 + vetor[1]**2 + vetor[2]**2)
		return (vetor[0]/sqt, vetor[1]/sqt, vetor[2]/sqt)

	def calForce(self, p1, p2):
		# 计算两个原子之间的势能，返回三个方向上的分量
		x = p1[0] - p2[0]
		y = p1[1] - p2[1]
		z = p1[2] - p2[2]		
		distance = numpy.sqrt(x**2 + y**2 + z**2)
		force = math.sin(min(distance,PI/2.0)**2)
		direct = self.direction((x,y,z))
		force_x = force * direct[0]
		force_y = force * direct[1]
		force_z = force * direct[2]

		return numpy.array((force_x, force_y, force_z))

	def calForceSum(self, p1, p_set):
		# 计算p1原子在所有其他原子影响下的势能总和，返回势能总和
		force = numpy.array((0,0,0))
		for p in p_set:
			force = force + self.calForce(p1,p)
		return force
		
	def calAccelerate(self, mole_num_min, mole_num_max):
		# 计算部分原子的加速度		
		#print '计算加速度(%d,%d)'%(mole_num_min, mole_num_max)		
		
		for j in range(mole_num_min, mole_num_max):
			p1 = (self.pos[j*DEMENTION], self.pos[j*DEMENTION+1], self.pos[j*DEMENTION+2])
			p_set = []
			for k in range(j):
				p_set.append((self.pos[k*DEMENTION], self.pos[k*DEMENTION+1], self.pos[k*DEMENTION+2]))
			for k in range(j+1,MOLECULAR_NUM):
				p_set.append((self.pos[k*DEMENTION], self.pos[k*DEMENTION+1], self.pos[k*DEMENTION+2]))
			# 返回合力在三个方向上的分力，分别除以原子质量得到分加速度，更新加速度
			a = self.calForceSum(p1,p_set)/M
			self.acc_dt[j*DEMENTION] = a[0]
			self.acc_dt[j*DEMENTION+1] = a[1]
			self.acc_dt[j*DEMENTION+2] = a[2]	

		if LOCK.acquire():
			global QUEUE, DATA_QUEUE
			QUEUE.put(1)
			DATA_QUEUE.put({mole_num_min/500: self.acc_dt[mole_num_min*DEMENTION:mole_num_max*DEMENTION]})
			LOCK.release()
	
	def updatePos(self, mole_num_min, mole_num_max):
		# 更新部分原子的位置
		#print '更新原子位置(%d,%d)'%(mole_num_min, mole_num_max)
		
		for j in range(mole_num_min, mole_num_max):
			self.pos[j*DEMENTION] = self.pos[j*DEMENTION] + self.vel[j*DEMENTION]*dt + 0.5*self.acc[j*DEMENTION]*dt**2
			self.pos[j*DEMENTION+1] = self.pos[j*DEMENTION+1] + self.vel[j*DEMENTION+1]*dt + 0.5*self.acc[j*DEMENTION+1]*dt**2
			self.pos[j*DEMENTION+2] = self.pos[j*DEMENTION+2] + self.vel[j*DEMENTION+2]*dt + 0.5*self.acc[j*DEMENTION+2]*dt**2
			for i in range(3):
				if self.pos[i+j*DEMENTION] > 10.0:
					self.pos[i+j*DEMENTION] = 10.0
				elif self.pos[i+j*DEMENTION] < 0.0:
					self.pos[i+j*DEMENTION] = 0.0

		if LOCK.acquire():
			global QUEUE, DATA_QUEUE
			QUEUE.put(1)			
			DATA_QUEUE.put({mole_num_min/500:self.pos[mole_num_min*DEMENTION:mole_num_max*DEMENTION]})
			LOCK.release()

	def updateVelocity(self, mole_num_min, mole_num_max):
		# 先更新部分原子的速度，再更新加速度
		#print '更新速度(%d,%d)'%(mole_num_min, mole_num_max)		
		
		for j in range(mole_num_min, mole_num_max):
			self.vel[j*DEMENTION] = self.vel[j*DEMENTION] + 0.5*(self.acc[j*DEMENTION] + self.acc_dt[j*DEMENTION])*dt
			self.vel[j*DEMENTION+1] = self.vel[j*DEMENTION+1] + 0.5*(self.acc[j*DEMENTION+1] + self.acc_dt[j*DEMENTION+1])*dt
			self.vel[j*DEMENTION+2] = self.vel[j*DEMENTION+2] + 0.5*(self.acc[j*DEMENTION+2] + self.acc_dt[j*DEMENTION+2])*dt		

		if LOCK.acquire():
			global QUEUE, DATA_QUEUE
			QUEUE.put(1)
			DATA_QUEUE.put({mole_num_min/500: self.vel[mole_num_min*DEMENTION:mole_num_max*DEMENTION]})
			LOCK.release()

	def updateState(self, spamwriter):
		global process_pool,lockCount
		while lockCount < 2:
			lockCount = lockCount + QUEUE.get()
			data = DATA_QUEUE.get()
			for n in data:
				self.vel[n*500*3:(n+1)*500*3] = data[n]
			time.sleep(0.2)
		lockCount = 0

		self.acc = self.acc_dt

		for process in process_pool:
			process.join()
		process_pool = []
		start_time = time.time()
		# 计算原子的加速度，更新到acc_dt	
		try:
			process_pool.append(Process(target=self.calAccelerate, args=(0,int(MOLECULAR_NUM*0.5))))
			process_pool.append(Process(target=self.calAccelerate, args=(int(MOLECULAR_NUM*0.5),MOLECULAR_NUM)))
			for process in process_pool:
				process.start()
		except Exception, e:
			print '加速度计算进程启动失败:',e

		while lockCount < 2:
			lockCount = lockCount + QUEUE.get()
			data = DATA_QUEUE.get()
			for n in data:
				self.acc_dt[n*500*3:(n+1)*500*3] = data[n]
			time.sleep(0.2)
		lockCount = 0
		print time.time()-start_time, 'secs to cal Accelerate'

		for process in process_pool:
			process.join()		
		process_pool = []
		# 更新原子的位置
		try:
			spamwriter.writerow(self.pos)
			process_pool.append(Process(target=self.updatePos, args=(0,int(MOLECULAR_NUM*0.5))))
			process_pool.append(Process(target=self.updatePos, args=(int(MOLECULAR_NUM*0.5),MOLECULAR_NUM)))
			for process in process_pool:
				process.start()
		except Exception, e:
			print '位置更新进程启动失败:',e

		while lockCount < 2:
			lockCount = lockCount + QUEUE.get()
			data = DATA_QUEUE.get()
			for n in data:
				self.pos[n*500*3:(n+1)*500*3] = data[n]
			time.sleep(0.2)
		lockCount = 0

		
		for process in process_pool:
			process.join()		
		process_pool = []
		# 更新原子的加速度和速度
		try:
			process_pool.append(Process(target=self.updateVelocity, args=(0,int(MOLECULAR_NUM*0.5))))
			process_pool.append(Process(target=self.updateVelocity, args=(int(MOLECULAR_NUM*0.5),MOLECULAR_NUM)))
			for process in process_pool:
				process.start()
		except Exception, e:
			print '速度更新进程启动失败:',e

if __name__ == '__main__':
	start_time = time.time()

	path = os.getcwd()	
	csvFile = open(os.path.join(path,'motionHistory.csv'),'w')
	spamwriter = csv.writer(csvFile,delimiter=',')

	msec = 0
	mms = MolecularMotionSimulator()
	while msec < 2000:		
		if MUTEX.acquire():
			print 'the',msec+1,'iteration'
			print 'atom 1 position (%f,%f,%f)'%(mms.pos[0],mms.pos[1],mms.pos[2])
			mms.updateState(spamwriter)
			msec = msec + 1
			MUTEX.release()
	csvFile.close()
	print time.time()-start_time, 'secs finished'


