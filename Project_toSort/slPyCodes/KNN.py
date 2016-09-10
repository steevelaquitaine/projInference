
# KNN.py
#
#     author: steeve laquitaine
#       date: 141225
#    purpose: K nearest neighbors classification
#description: machine learning

#load package
import csv
import LoadData
import euclideanDistance
import getNeighbors

#load data
with open('iris.data', 'rb') as csvfile:
	
	# read data
	lines = csv.reader(csvfile)
	for row in lines:
		print ', '.join(row)


# split data fosaar training KNN and test KNN
trainingSet=[]
testSet    =[]
LoadData.LoadData('iris.data', 0.66, trainingSet, testSet)
print 'Train: ' + repr(len(trainingSet))
print 'Test:  ' + repr(len(testSet))


# Similarity
data1 = [2, 2, 2, 'a']
data2 = [4, 4, 4, 'b']
distance = euclideanDistance.euclideanDistance(data1, data2, 3)
print 'Distance: ' + repr(distance)


# Neighborsfire
trainSet = [[2, 2, 2, 'a'], [4, 4, 4, 'b']]
testInstance = [5, 5, 5]
k = 4
neighbors = getNeighbors.getNeighbors(trainSet, testInstance, 1)
print(neighbors)

import operator
def getResponse(neighbors):

	# create dictionary {response: vote}
	classVotes = {}

	# each neighbor
	for x in range(len(neighbors)):

		# neighbor response (last column)
		response = neighbors[x][-1]

		# add neighbor's vote to possible responses 
		if response in classVotes:
			classVotes[response] += 1
		else:
			classVotes[response] = 1
			
    # most voted first
	sortedVotes = sorted(classVotes.iteritems(), key=operator.itemgetter(1), reverse=True)
	return sortedVotes[0][0]


neighbors = [[1,1,1,'a'], [2,2,2,'a'], [3,3,3,'b']]
response = getResponse(neighbors)
print(response)
