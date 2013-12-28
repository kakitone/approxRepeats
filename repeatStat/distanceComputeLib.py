from operator import itemgetter
import matplotlib.pyplot as plt
import math
import numpy as np
import os
import csv
import bisect
import random

def hammingDistance(str1, str2, k):
    count  =0
    for index in range(k):
        if str1[index] != str2[index]:
            count = count +1
    return count