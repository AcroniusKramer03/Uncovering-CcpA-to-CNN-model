'''
Author: Anne de Jong 
Adapted script from pppv2_CNN_Model_Hung_Anne.py to make it suitable for building TF/promoter models
Date: 17/04/2024

python3 ProPr_makeCNNmodel_from_SeqList.py -sessiondir examples 

'''


#import os
import sys
import pandas as pd
import numpy as np
import argparse
import re
import random
import tensorflow as tf
import time

from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.model_selection import train_test_split
from tensorflow.keras.layers import Conv1D, Dense, MaxPooling1D, Flatten, Bidirectional
from tensorflow.keras.models import Sequential
from tensorflow.keras.callbacks import ModelCheckpoint, History, ReduceLROnPlateau, EarlyStopping
from tensorflow.keras.models import load_model
from tensorflow.python.keras.utils.data_utils import Sequence
from sklearn.metrics import confusion_matrix
from matplotlib import pyplot as plt

# ProPr libs
from models import *


random.seed(0)

# ---------------------------------------------------------------- parse parameters ---------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Train the model')
parser.add_argument('-sessiondir', dest='sessiondir', help='Session Dir [.]', nargs='?', default='.')
parser.add_argument('-tr', dest='training', help='training_seqs filename in sessiondir [training_seqs.txt]', default='training_seqs.txt')
parser.add_argument('-bg', dest='background', help='background_seqs filename in sessiondir [background_seqs.txt]', default='background_seqs.txt')
parser.add_argument('-dataprep', dest='dataPrep', help='Whether to perform data preparation or not [true]', default='true')  # Hung
parser.add_argument('-training', dest='training', help='Whether to perform model training and validation [true]', default='true')  # Hung
parser.add_argument('-modelname',dest='modelName', help='Name of the model without extension [cnn_lstm]', nargs='?', default='cnn_lstm')  # Hung
parser.add_argument('-testset',dest='testSet', help='Percent of the data for the test set [0.1]', nargs='?', default=0.1) 
parser.add_argument('-epochs',dest='epochs', help='Number of training epochs [11]', nargs='?', default=11)  # Hung
parser.add_argument('-batchsize',dest='batchSize', help='Batch size 32,64,128 [64]', nargs='?', default=64)  # Hung
parser.add_argument('-train', dest='TrainOutfile', help='Prefix of theFrank file [Train]', nargs='?', default='Train')
parser.add_argument('-frank', dest='FrankFile', help='Output filename prefix for Frank; prefix.fna and prefix.gff [Frank]', nargs='?', default='Frank')
parser.add_argument('--version', action='version', version='Anne de Jong, version 1.0, Feb 2020')

args = parser.parse_args()

# Be sure the have real numbers
args.epochs = int(args.epochs)
args.batchSize = int(args.batchSize)
args.testSet = float(args.testSet)

gff_header = ["chrom","db", "type", "start", "end", "name", "strand", "score", "description"]
bed_header = ["chrom","start", "end", 'name','score', "strand"]

''' --------------------------------------  DEFINE FUNCTIONS ---------------------------------------------------'''

def Anne_one_hot_encode(seq):
  mapping = dict(zip("ACGT", range(4)))
  seq2 = [mapping[i] for i in seq]
  return np.eye(4)[seq2]

# Hung (2): for printing memory usage of all variables
def sizeof_fmt(num, suffix='B'):
  ''' by Fred Cirera,  https://stackoverflow.com/a/1094933/1870254, modified'''
  for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
    if abs(num) < 1024.0:
      return "%3.1f %s%s" % (num, unit, suffix)
    num /= 1024.0
  return "%.1f %s%s" % (num, 'Yi', suffix)


def write_log(S):
	# write to console and logfile
	print(S)
	f = open(args.sessiondir+'/ProPr_makeCNNmodel.log', 'a')
	f.write(S + '\n')
	f.close()


''' ======================================   MAIN ================================================'''
 

if args.dataPrep == 'true':
	''' ====  Load training sets ====='''
	with open(args.sessiondir +'/'+args.training, 'r')   as f: 
		training_list   = [line.strip() for line in f.readlines()]
	with open(args.sessiondir +'/'+args.background, 'r') as f: 
		background_list = [line.strip() for line in f.readlines()]
	
	write_log('Number of training seqs:   ' + str(len(training_list)))
	write_log('Number of background seqs: ' + str(len(background_list)))
	
	''' ======================================  TRAIN THE MODEL ================================================'''
	write_log('Training the model...')
	
	# 1. Seperate out a test set
	write_log('Seperate out a test set...')
	test_set_fraction = 0.1
	test_set_percentage = 100 / (100 * test_set_fraction)
	write_log('Using ' +str(test_set_percentage) +' percent of the data for the test set\n')
	
	training_sequences = []
	training_response = []
	test_sequences = []
	test_response = []
		
	# get test sequences based on percentage of test_set_percentage  and add the '1'  response value
	for i in range (0,len(training_list)):
		if i % test_set_percentage == 0:
			test_sequences.append(training_list[i])
			test_response.append(1)
		else:
			training_sequences.append(training_list[i])
			training_response.append(1)
	
	# get the background fraction and add the '0  response value
	for i in range (0,len(background_list)):
		if i % test_set_percentage == 0:
			test_sequences.append(background_list[i])
			test_response.append(0)
		else:
			training_sequences.append(background_list[i])
			training_response.append(0)
	
	write_log('Number of sequences in training set:  '+str(len(training_sequences)))
	write_log('Number of sequences in test set:      '+str(len(test_sequences)))
	write_log('===========')
	
	
	# 2. Get sequences ready for training as features
	write_log('Encode sequences to an image...')
	
	# The LabelEncoder encodes a sequence of bases as a sequence of integers.
	# OLD integer_encoder = LabelEncoder()
	
	# The OneHotEncoder converts an array of integers to a sparse matrix where
	one_hot_encoder = OneHotEncoder(categories='auto')
	
	write_log('Convert training sequences to image with one_hot_encode')
	train_features = []
	for sequence in training_sequences: train_features.append(Anne_one_hot_encode(sequence))
	np.set_printoptions(threshold=40)
	train_features = np.stack(train_features)
	
	write_log('Convert test sequences to image with one_hot_encode')
	test_features = []
	for sequence in test_sequences:	test_features.append(Anne_one_hot_encode(sequence))
	np.set_printoptions(threshold=40)
	test_features = np.stack(test_features)
	
	
	# 3. Get Responses ready for training as labels
	train_labels = training_response
	one_hot_encoder = OneHotEncoder(categories='auto')
	train_labels = np.array(train_labels).reshape(-1, 1)
	train_labels = one_hot_encoder.fit_transform(train_labels).toarray()
	test_labels = test_response
	one_hot_encoder = OneHotEncoder(categories='auto')
	test_labels = np.array(test_labels).reshape(-1, 1)
	test_labels = one_hot_encoder.fit_transform(test_labels).toarray()
	
	# 4. OAnne to save disk space Save the numpy data for testing purposes
	np.save(args.sessiondir + '/train_features.npy', train_features)
	np.save(args.sessiondir + '/train_labels.npy' , train_labels)
	np.save(args.sessiondir + '/test_features.npy', test_features)
	np.save(args.sessiondir + '/test_labels.npy', test_labels)

	# Remove test_features and test_labels to clear up memory
	del test_features
	del test_labels
	
else:
	# For testing load previous numpy data to skip this part
	train_features = np.load(args.sessiondir + '/train_features.npy') 
	train_labels   = np.load(args.sessiondir + '/train_labels.npy') 
	test_features  = np.load(args.sessiondir + '/test_features.npy') 
	test_labels    = np.load(args.sessiondir + '/test_labels.npy') 


''' ============ model architecture ==================='''

# 5. selecting the model architecture
input_shape =(train_features.shape[1], 4)
print("args.modelName:", args.modelName)
model = create_model(input_shape=input_shape, model_name=args.modelName)
model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['binary_accuracy'])
print(model.summary())
model_filepath = args.sessiondir + '/' + args.modelName + '.keras'
write_log('Model will be saved to: ' + model_filepath)

# Hung 
checkpoint = ModelCheckpoint(filepath=model_filepath, save_weights_only=False, monitor='val_loss', mode='min', save_best_only=True)
reduce_lr = ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=5, min_lr=0.00001)
early_stopping = EarlyStopping(monitor='val_loss', patience=8, mode='min', restore_best_weights=True)

# Train_validation split
if 'embed' in args.modelName:
	train_features = np.argmax(train_features, axis=-1)  # map (N, 71, 4) to (N, 71)
X_train, X_val, y_train, y_val = train_test_split(train_features, train_labels, test_size=0.25, random_state=0)
print("X_train.shape:", X_train.shape)
print("X_val.shape:"  , X_val.shape)
print("y_train.shape:", y_train.shape)
print("y_val.shape:"  , y_val.shape)

# Print memory usage of all variables
for name, size in sorted(((name, sys.getsizeof(value)) for name, value in locals().items()),key= lambda x: -x[1])[:10]:
	print("{:>30}: {:>8}".format(name, sizeof_fmt(size)))

history = model.fit(train_features, train_labels, epochs=args.epochs, verbose=2, validation_split=0.25)





''' ============ save model and report performance ==================='''

history_df = pd.DataFrame(history.history)
history_df.to_csv(args.sessiondir+'/'+args.modelName+'.history.txt', index = True, sep ='\t')

plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'val'], loc='upper left')
#plt.show()
plt.savefig(args.sessiondir+'/'+'model_loss_plot.png')

# 6. Save / Load the model
write_log('Model saved to: ' + model_filepath) 
model.save(model_filepath)

# 7. Validation using test set
def Validation_report(features, labels, name):
	# F1 = A measure that combines precision and recall is the harmonic mean of precision and recall, the traditional F-measure or balanced F-score:
	# F1 = 2 * (precision * recall) / (precision + recall)
	# precision = TP / TP + FP
	# recall = TP / TP + FN
	predicted_labels = model.predict(features)
	cm = confusion_matrix(np.argmax(labels, axis=1), np.argmax(predicted_labels, axis=1))
	TP = cm[1][1]
	FP = cm[0][1]
	FN = cm[1][0]
	precision = TP / ( TP + FP )
	recall = TP / (TP + FN)
	F1 = 2 * (precision * recall) / (precision + recall)

	write_log('Confusion matrix of '+ name+'\n' + str(cm) + '| '+args.sessiondir)
	write_log('Precision  TP / TP + FP                              '+str(round(precision,2)))
	write_log('Recall     TP / TP + FN                              '+str(round(recall,2)))
	write_log('F1 = 2 * (precision * recall) / (precision + recall) '+str(round(F1,2)))

	# Confusion matrix info is in the log, but one extra clean file:
	f = open(args.sessiondir+'/'+"Confusion_matrix.txt", "a")
	f.write(args.sessiondir+ '\n' )
	f.write('Confusion matrix of '+ name+'\n' + str(cm) + '\n' )
	f.write('Precision  TP / TP + FP                              '+str(round(precision,2))+ '\n')
	f.write('Recall     TP / TP + FN                              '+str(round(recall,2))+ '\n')
	f.write('F1 = 2 * (precision * recall) / (precision + recall) '+str(round(F1,2))+ '\n\n')
	f.write('Precision\tRecall\tF1\n')
	f.write(str(round(precision,2))+'\t'+str(round(recall,2))+'\t'+str(round(F1,2))+'\n\n')
	f.close() 
	
random.seed(0)
Validation_report(train_features, train_labels, 'Training_set') 

# 8. Done
write_log("============\n|   DONE   |\n============\n")








''' ================================ OUTPUT TEST PREDICTIONS FOR FURTHER INVESTIGATION =========================================='''

#    correct_predictions = []
#    false_predictions = []
#    
#    for i in range(0, len(test_features)):
#        if test_labels[i][1] == 1 and round(predicted_labels[i][1]) == 1:
#            correct_predictions.append(test_features[i])
#        if test_labels[i][1] == 1 and round(predicted_labels[i][1]) == 0:
#            false_predictions.append(test_features[i])
#    
#    ID_correct = []
#    value = 1
#    for i in correct_predictions:
#        if i in promoter_list:
#            header = promoter_list.index(i) - 1
#            ID_correct.append(promoter_list[header])
#        else:
#            ID_correct.append('>NONPROMOTER|'+str(value))
#            value +=1
#    
#    
#    ID_false = []
#    value = 1
#    for i in false_predictions:
#        if i in promoter_list:
#            header = promoter_list.index(i) - 1
#            ID_false.append(promoter_list[header])
#        else:
#            ID_false.append('>NONPROMOTER|'+str(value))
#            value +=1
#    
#    
#    f = open(args.sessiondir+'/correct_predictions.txt', 'w')
#    for i in range(0, len(correct_predictions)):
#        f.write(str(ID_correct[i]))
#        f.write('\n')
#        f.write(str(correct_predictions[i]))
#        f.write('\n')
#    f.close()
#    
#    f =   open(args.sessiondir+'/false_predictions.txt', 'w')
#    for i in range(0, len(false_predictions)):
#        f.write(str(ID_false[i]))
#        f.write('\n')
#        f.write(str(false_predictions[i]))
#        f.write('\n')
#    f.close()
#    

# cm = cm.astype('float') / cm.sum(axis = 1)[:, np.newaxis]
#
# plt.imshow(cm, cmap=plt.cm.Blues)
# plt.title('Normalized confusion matrix')
# plt.colorbar()
# plt.xlabel('Predicted label')
# plt.ylabel('True label')
# plt.xticks([0, 1]); plt.yticks([0, 1])
# plt.grid('off')
# for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
#     plt.text(j, i, format(cm[i, j], '.2f'),
#              horizontalalignment='center',
#              color='white' if cm[i, j] > 0.5 else 'black')
# plt.show()
