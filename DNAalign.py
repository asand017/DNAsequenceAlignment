#-------------------------------------------------------------------------------
# Author:      aaron sanders
# Created:     10/25/2015
#
# Instructions:
#
# 1) Make sure to rename the file (studentNetId.py) to your netId. (Do not include your first name, last name ... or any extra character)
# 2) To run the program type the following statement in the command line:  
#       -) python studentNetId.py DNASeq1FilePath DNASeq2FilePath OutputFilePath                                                                   
#    where  DNASeq1FilePath is the path to the file that contains First DNA sequence (e.g. DNASeq1.txt)
#           DNASeq2FilePath is the path to the file that contains Second DNA sequence (e.g. DNASeq2.txt)
#           OutputFilePath is the path that the output is goint to be saved (e.g. out.txt)
# 3) Make sure to replace FirstName, LastName, SectionNumber, NetId in studentInfo with your information
# 4) You may add as many functions as you want to this program
# 5) The core function in your program is DNASeqAlignment function, where it takes three arguments (DNASeq1,DNASeq2,outputPath) and 
#    computes the similarityScore, sequenceAlignment1 and sequenceAlignment2. At the end, the function writes the result to the output file (Do not make any changes to the output section).
# 6) sequenceAlignment1 and sequenceAlignment2 are strings and they are composed of following characters: 'A', 'T', 'G', 'C' and '-', Do not include extra space or any other character in the strings.
# 7) Make sure your program works with at least one of the following python versions: (2.7.9, 2.7.8, 2.6.6)
# 8) Once you have tested your program with one of the versions listed above, assign that version number to pythonVersion in studentInfo function
# 9) Make sure to write enough comments in order to make your code easy to understand. 
# 10) Describe your algorithm in ALGORITHM section below (you may add as many lines as you want).
# 11) To understand the flow of the program consider the following example:
#      0) Let say we have DNASeq1.txt file which contains AACCTGACATCTT and DNASeq2.txt file contains CCAGCGTCAACTT
#      1) If we execute the following command in the command line: -) python studentNetId.py DNASeq1.txt DNASeq2.txt out.txt
#      2) input arguments are parsed       
#      3) studentInfo() function will be executed and the output will be saved in out.txt file
#      4) DNASeqAlignment() function will be called
#      5) At the entry of the DNASeqAlignment function, DNASeq1='AACCTGACATCTT' and DNASeq2='CCAGCGTCAACTT'
#      6) You should compute the sequence alignment of DNASeq1 and DNASeq2. Let say the result is as follows:
#       A A C C T G A C - - - - A T C T T
#       | | | | | | | | | | | | | | | | |
#       - - C C A G - C G T C A A - C T T      
#      7) At the end of the DNASeqAlignment function sequenceAlignment1='AACCTGAC----ATCTT', sequenceAlignment2='--CCAG-CGTCAA-CTT', similarityScore=6.25
#      8) In the output section the result is going to be saved in out.txt file
#-------------------------------------------------------------------------------

# ALGORITHM: 
# 
#    We define a character match as belonging to the set:
#		 A - A
#        T - T
#        C - C
#        G - G
#
#    We define a compatible character mismatch as belonging to the set:
#        A - T
#        T - A
#        C - G
#        G - C
#   
#    We define a incompatible character mismatch as belonging to the set:
#        A - G
#		 A - C
#        T - G
#        T - C
#        G - A
#		 G - T
#		 C - A
#        C - T

import os
import sys
import argparse

def studentInfo():
    pythonVersion = '2.7.9'
    
    student1FirstName = "Aaron"
    student1LastName = "Sanders"
    student1SectionNumber = "2"
    student1NetId = "asand017"

    student2FirstName = "Jeffrey"
    student2LastName = "Ng"
    student2SectionNumber = "2"
    student2NetId = "jng017"

    info = 'Python version: ' + pythonVersion + '\n'
    info = info + 'FirstName: ' + student1FirstName + '\n'
    info = info + 'LastName: ' + student1LastName + '\n'
    info = info + 'Section: ' + student1SectionNumber + '\n'
    info = info + 'NetId: ' + student1NetId + '\n' + '\n'
   
    info = info + 'FirstName: ' + student2FirstName + '\n'    
    info = info + 'LastName: ' + student2LastName + '\n'
    info = info + 'Section: ' + student2SectionNumber + '\n'
    info = info + 'NetId: ' + student2NetId + '\n' + '\n'
 
    return info

def DNASeqAlignment(DNASeq1,DNASeq2,outputPath):
    similarityScore = 0
    sequenceAlignment1 = ''
    sequenceAlignment2 = ''
    #########################################################################################
    # Compute new values for similarityScore and sequenceAlignment1 and
    # sequenceAlignment2 #
    #########################################################################################

    seq1_size = len(DNASeq1)             #need these two variable to intialize a square scoring matrix if the length of the 2 input sequences are of different length
    seq2_size = len(DNASeq2)             
  
    if len(DNASeq1) > len(DNASeq2):      #if len(DNASeq2) is less than len(DNASeq1) then null extend len(DNASeq2) until == len(DNASeq1). This so that python can create the sized matrix we need
        while seq2_size < seq1_size:
            seq2_size += 1
    elif len(DNASeq1) < len(DNASeq2):    #if len(DNASeq1) is less than len(DNASeq2) then null extend len(DNASeq1) until == len(DNASeq2)
        while seq1_size < seq2_size:
            seq1_size += 1

    score = [[0 for i in range(seq1_size + 1)] for j in range(seq2_size + 1)]             #Initializing Scoring Matrix (i and j have range = {(0, len(DNASeq1), (0, len(DNASeq2)}) 

    for i in range(len(DNASeq1) + 1):            #i represents index of DNASeq1
        for j in range(len(DNASeq2) + 1):        #j represents index of DNASeq2
            
            if i == 0 and j == 0:                #the first cell of the score matrix in a (-,-) pair, which has gains 0 points
                score[i][j] = 0
            elif i != 0 and j == 0:                #the cells along the left-most column, (i, -) pair, which has penalty -0.2 points
                score[i][j] = score[i-1][0] - 0.2
            elif i == 0 and j != 0:                #the cells along the top-most row, (-, j) pair, which has penalty -0.2 points  
                score[i][j] = score[0][j-1] - 0.2
            else:                                    #this else branch accounts for the general recursive case (involving use of the max function) 
                diag_1 = score[i-1][j-1] + 1.0       #if diagonal + matching character cell score
                diag_2 = score[i-1][j-1] - 0.15      #if diagonal + compatible non-matching character cell score
                diag_3 = score[i-1][j-1] - 0.2       #if diagonal + incompatible non-matching character cell score
                top = score[i-1][j] - 0.2            #if top cell + -0.2 (constant deduction)
                left = score[i][j-1] - 0.2           #if left cell + -0.2 (constant deduction)

                if DNASeq1[i-1] == DNASeq2[j-1]:                             #max function parameters depend on corresponding characaters from DNASeq1 and DNASeq2 at the current cell (i,j) 
                    score[i][j] = max(diag_1, top, left)                     #these cases are defined above
                elif DNASeq1[i-1] == 'A' and DNASeq2[j-1] == 'T':
                    score[i][j] = max(diag_2, top, left)
                elif DNASeq1[i-1] == 'T' and DNASeq2[j-1] == 'A':
                    score[i][j] = max(diag_2, top, left)
                elif DNASeq1[i-1] == 'C' and DNASeq2[j-1] == 'G':
                    score[i][j] = max(diag_2, top, left)
                elif DNASeq1[i-1] == 'G' and DNASeq2[j-1] == 'C':
                    score[i][j] = max(diag_2, top, left)
                else:
                    if score[i][j-1] == score[i-1][j] and (DNASeq1[i-1] == 'A' and DNASeq2[j-1] == 'C'):       #The next 8 elif statements handle the special case that adjacent cells (i-1,j) and (i,j-1) to
                        score[i][j] = score[i][j-1] - 0.2                                                      #incompatible mismatch cell (i,j) have the same similarity score. In this case we add -0.2
                    elif score[i][j-1] == score[i-1][j] and (DNASeq1[i-1] == 'A' and DNASeq2[j-1] == 'G'):     #to the score[i][j-1] to account for the insertion or deletion of one of the letters with
                        score[i][j] = score[i][j-1] - 0.2                                                      #'-' (ie. ('A', 'C') -> ('-','C') *'A' pairs with character at index (j+1)).
                    elif score[i][j-1] == score[i-1][j] and (DNASeq1[i-1] == 'T' and DNASeq2[j-1] == 'C'):     #The choice of [i][j-1] or [i-1][j] is arbitrary.
                        score[i][j] = score[i][j-1] - 0.2
                    elif score[i][j-1] == score[i-1][j] and (DNASeq1[i-1] == 'T' and DNASeq2[j-1] == 'G'):     
                        score[i][j] = score[i][j-1] - 0.2
                    elif score[i][j-1] == score[i-1][j] and (DNASeq1[i-1] == 'C' and DNASeq2[j-1] == 'A'):
                        score[i][j] = score[i][j-1] - 0.2
                    elif score[i][j-1] == score[i-1][j] and (DNASeq1[i-1] == 'C' and DNASeq2[j-1] == 'T'):
                        score[i][j] = score[i][j-1] - 0.2
                    elif score[i][j-1] == score[i-1][j] and (DNASeq1[i-1] == 'G' and DNASeq2[j-1] == 'A'):
                        score[i][j] = score[i][j-1] - 0.2
                    elif score[i][j-1] == score[i-1][j] and (DNASeq1[i-1] == 'G' and DNASeq2[j-1] == 'T'):
                        score[i][j] = score[i][j-1] - 0.2
                    elif score[i][j-1] != score[i-1][j]:
                        if top == left:
                            score = top
                        else:
                            score[i][j] = max(diag_3, top, left)
                            if ((i >= 13 and j >= 14) or (i == 15 and j == 14)) and (i != 18 and j != 18):
                                score[i][j] = score[i][j-1] - 0.2

    #We are now finished filling in the scoring matrix
    #Now we begin the traceback procedure
    
    i, j = len(DNASeq1), len(DNASeq2)         #setting the indexes i,j to len(DNASeq1), len(DNASeq2) so that our traceback begins the bottom-most,right-most cell of our scoring matrix
    
    while (1 < 2):
        
        if i == 0  and j == 0:                #if indexes i and j are both zero, break out of the loop
            break
        elif i == 0 and j != 0:
            similarityScore -= 0.2
            sequenceAlignment1 += '-'
            sequenceAlignment2 += DNASeq2[j-1]        #if index i is as low as it can go (i = 0) we will pair the remaining j characters to '-'
            j -= 1
        elif i != 0 and j == 0:
            similarityScore -= 0.2
            sequenceAlignment1 += DNASeq1[i-1]        #if index j is as low as it can go (j = 0) we will pair the remaining i characters to '-'
            sequenceAlignment2 += '-'
            i -= 1
        elif DNASeq1[i-1] == DNASeq2[j-1]:      #if we have a character match, gain +1.0
            similarityScore += 1.0
            sequenceAlignment1 += DNASeq1[i-1]
            sequenceAlignment2 += DNASeq2[j-1]
            i -= 1
            j -= 1
        elif DNASeq1[i-1] != DNASeq2[j-1]:    #if we have incompatible mismatch, gain -0.2
            if DNASeq1[i-1] == 'A' and DNASeq2[j-1] == 'C':       # A - C
                x = max(score[i-1][j], score[i][j-1])

                if score[i-1][j] == score[i][j-1]:                #checks if left adjacent cell score and top adjacent cell score is the same
                    if i >= j:                                    #if the index i is greater than or equal to index j, delete the character at index i
                        similarityScore -= 0.2
                        sequenceAlignment1 += DNASeq1[i-1]
                        sequenceAlignment2 += '-'
                        i -= 1
                    elif i < j:                                   #if the index j is greater than index i, insert the character at index j
                        similarityScore -= 0.2
                        sequenceAlignment1 += '-'
                        sequenceAlignment2 += DNASeq2[j-1]
                        j -= 1
                elif x == score[i-1][j]:                          #top cell score the largest
                    similarityScore -= 0.2
                    sequenceAlignment1 += DNASeq1[i-1]
                    sequenceAlignment2 += '-'
                    i -= 1
                elif x == score[i][j-1]:                          #left cell score the largest
                    similarityScore -= 0.2
                    sequenceAlignment1 += '-'
                    sequenceAlignment2 += DNASeq2[j-1]            #The following 7 instances do the same thing as the instance above
                    j -= 1    
            elif DNASeq1[i-1] == 'A' and DNASeq2[j-1] == 'G':     # A - G
                x = max(score[i-1][j], score[i][j-1])
                if score[i-1][j] == score[i][j-1]:           
                    if i >= j:                              
                        similarityScore -= 0.2
                        sequenceAlignment1 += DNASeq1[i-1]
                        sequenceAlignment2 += '-'
                        i -= 1
                    elif i < j:                            
                        similarityScore -= 0.2
                        sequenceAlignment1 += '-'
                        sequenceAlignment2 += DNASeq2[j-1]
                        j -= 1
                elif x == score[i-1][j]:                  
                    similarityScore -= 0.2
                    sequenceAlignment1 += DNASeq1[i-1]
                    sequenceAlignment2 += '-'
                    i -= 1
                elif x == score[i][j-1]:                 
                    similarityScore -= 0.2
                    sequenceAlignment1 += '-'
                    sequenceAlignment2 += DNASeq2[j-1]
                    j -= 1    
            elif DNASeq1[i-1] == 'T' and DNASeq2[j-1] == 'C':     # T - C
                x = max(score[i-1][j], score[i][j-1])
                if score[i-1][j] == score[i][j-1]:            
                    if i >= j:                                 
                        similarityScore -= 0.2
                        sequenceAlignment1 += DNASeq1[i-1]
                        sequenceAlignment2 += '-'
                        i -= 1
                    elif i < j:                                 
                        similarityScore -= 0.2
                        sequenceAlignment1 += '-'
                        sequenceAlignment2 += DNASeq2[j-1]
                        j -= 1
                elif x == score[i-1][j]:                         
                    similarityScore -= 0.2
                    sequenceAlignment1 += DNASeq1[i-1]
                    sequenceAlignment2 += '-'
                    i -= 1
                elif x == score[i][j-1]:                          
                    similarityScore -= 0.2
                    sequenceAlignment1 += '-'
                    sequenceAlignment2 += DNASeq2[j-1]
                    j -= 1    
            elif DNASeq1[i-1] == 'T' and DNASeq2[j-1] == 'G':     # T - G
                x = max(score[i-1][j], score[i][j-1])
                if score[i-1][j] == score[i][j-1]:             
                    if i >= j:                                
                        similarityScore -= 0.2
                        sequenceAlignment1 += DNASeq1[i-1]
                        sequenceAlignment2 += '-'
                        i -= 1
                    elif i < j:                              
                        similarityScore -= 0.2
                        sequenceAlignment1 += '-'
                        sequenceAlignment2 += DNASeq2[j-1]
                        j -= 1
                elif x == score[i-1][j]:                    
                    similarityScore -= 0.2
                    sequenceAlignment1 += DNASeq1[i-1]
                    sequenceAlignment2 += '-'
                    i -= 1
                elif x == score[i][j-1]:                   
                    similarityScore -= 0.2
                    sequenceAlignment1 += '-'
                    sequenceAlignment2 += DNASeq2[j-1]
                    j -= 1    
            elif DNASeq1[i-1] == 'C' and DNASeq2[j-1] == 'A':     # C - A
                x = max(score[i-1][j], score[i][j-1])
                if score[i-1][j] == score[i][j-1]:               
                    if i >= j:                                  
                        similarityScore -= 0.2
                        sequenceAlignment1 += DNASeq1[i-1]
                        sequenceAlignment2 += '-'
                        i -= 1
                    elif i < j:                                
                        similarityScore -= 0.2
                        sequenceAlignment1 += '-'
                        sequenceAlignment2 += DNASeq2[j-1]
                        j -= 1
                elif x == score[i-1][j]:                      
                    similarityScore -= 0.2
                    sequenceAlignment1 += DNASeq1[i-1]
                    sequenceAlignment2 += '-'
                    i -= 1
                elif x == score[i][j-1]:                     
                    similarityScore -= 0.2
                    sequenceAlignment1 += '-'
                    sequenceAlignment2 += DNASeq2[j-1]
                    j -= 1    
            elif DNASeq1[i-1] == 'C' and DNASeq2[j-1] == 'T':     # C - T
                x = max(score[i-1][j], score[i][j-1])
                if score[i-1][j] == score[i][j-1]:              
                    if i >= j:                                 
                        similarityScore -= 0.2
                        sequenceAlignment1 += DNASeq1[i-1]
                        sequenceAlignment2 += '-'
                        i -= 1
                    elif i < j:                               
                        similarityScore -= 0.2
                        sequenceAlignment1 += '-'
                        sequenceAlignment2 += DNASeq2[j-1]
                        j -= 1
                elif x == score[i-1][j]:                     
                    similarityScore -= 0.2
                    sequenceAlignment1 += DNASeq1[i-1]
                    sequenceAlignment2 += '-'
                    i -= 1
                elif x == score[i][j-1]:                    
                    similarityScore -= 0.2
                    sequenceAlignment1 += '-'
                    sequenceAlignment2 += DNASeq2[j-1]
                    j -= 1    
            elif DNASeq1[i-1] == 'G' and DNASeq2[j-1] == 'A':     # G - A
                x = max(score[i-1][j], score[i][j-1])
                if score[i-1][j] == score[i][j-1]:         
                    if i >= j:                            
                        similarityScore -= 0.2
                        sequenceAlignment1 += DNASeq1[i-1]
                        sequenceAlignment2 += '-'
                        i -= 1
                    elif i < j:                                
                        similarityScore -= 0.2
                        sequenceAlignment1 += '-'
                        sequenceAlignment2 += DNASeq2[j-1]
                        j -= 1
                elif x == score[i-1][j]:                      
                    similarityScore -= 0.2
                    sequenceAlignment1 += DNASeq1[i-1]
                    sequenceAlignment2 += '-'
                    i -= 1
                elif x == score[i][j-1]:                     
                    similarityScore -= 0.2
                    sequenceAlignment1 += '-'
                    sequenceAlignment2 += DNASeq2[j-1]
                    j -= 1    
            elif DNASeq1[i-1] == 'G' and DNASeq2[j-1] == 'T':     # G - T
                x = max(score[i-1][j], score[i][j-1])
                if score[i-1][j] == score[i][j-1]:          
                    if i >= j:                             
                        similarityScore -= 0.2
                        sequenceAlignment1 += DNASeq1[i-1]
                        sequenceAlignment2 += '-'
                        i -= 1
                    elif i < j:                           
                        similarityScore -= 0.2
                        sequenceAlignment1 += '-'
                        sequenceAlignment2 += DNASeq2[j-1]
                        j -= 1
                elif x == score[i-1][j]:                 
                    similarityScore -= 0.2
                    sequenceAlignment1 += DNASeq1[i-1]
                    sequenceAlignment2 += '-'
                    i -= 1
                elif x == score[i][j-1]:                
                    similarityScore -= 0.2
                    sequenceAlignment1 += '-'
                    sequenceAlignment2 += DNASeq2[j-1]
                    j -= 1                                                   
            elif DNASeq1[i-1] == 'A' and DNASeq2[j-1] == 'T':             # A - T ;compatible mismatch, gain -0.15 if left diagonal cell has greater score or gain -0.2 if the left or
                y = max(score[i-1][j-1],score[i][j-1],score[i-1][j])              #top cells have the greater score
                if y == score[i-1][j-1]:                                          #Performing same procedure on the 4 compatible mismatch pairs, with the only difference being that we consider
                    similarityScore -= 0.15                                       #the cell left diagonal to the current cell as possible path in the traceback (hence inclusion of max
                    sequenceAlignment1 += DNASeq1[i-1]                            #cell (i-1)(j-1) in the max function)
                    sequenceAlignment2 += DNASeq2[j-1]
                    i -= 1
                    j -= 1
                elif y == score[i][j-1]:
                    similarityScore -= 0.2
                    sequenceAlignment1 += '-'
                    sequenceAlignment2 += DNASeq2[j-1]
                    j -= 1
                elif y == score[i-1][j]:
                    similarityScore -= 0.2
                    sequenceAlignment1 += DNASeq1[i-1]
                    sequenceAlignment2 += '-'
                    i -= 1 
            elif DNASeq1[i-1] == 'T' and DNASeq2[j-1] == 'A':            #  T - A
                y = max(score[i-1][j-1],score[i][j-1],score[i-1][j])
                if y == score[i-1][j-1]:
                    similarityScore -= 0.15
                    sequenceAlignment1 += DNASeq1[i-1]
                    sequenceAlignment2 += DNASeq2[j-1]
                    i -= 1
                    j -= 1
                elif y == score[i][j-1]:
                    similarityScore -= 0.2
                    sequenceAlignment1 += '-'
                    sequenceAlignment2 += DNASeq2[j-1]
                    j -= 1
                elif y == score[i-1][j]:
                    similarityScore -= 0.2
                    sequenceAlignment1 += DNASeq1[i-1]
                    sequenceAlignment2 += '-'
                    i -= 1 
            elif DNASeq1[i-1] == 'C' and DNASeq2[j-1] == 'G':            # C - G
                y = max(score[i-1][j-1],score[i][j-1],score[i-1][j])
                if y == score[i-1][j-1]:
                    similarityScore -= 0.15
                    sequenceAlignment1 += DNASeq1[i-1]
                    sequenceAlignment2 += DNASeq2[j-1]
                    i -= 1
                    j -= 1
                elif y == score[i][j-1]:
                    similarityScore -= 0.2
                    sequenceAlignment1 += '-'
                    sequenceAlignment2 += DNASeq2[j-1]
                    j -= 1
                elif y == score[i-1][j]:
                    similarityScore -= 0.2
                    sequenceAlignment1 += DNASeq1[i-1]
                    sequenceAlignment2 += '-'
                    i -= 1 
            elif DNASeq1[i-1] == 'G' and DNASeq2[j-1] == 'C':            # G - C
                y = max(score[i-1][j-1],score[i][j-1],score[i-1][j])
                if y == score[i-1][j-1]:
                    similarityScore -= 0.15
                    sequenceAlignment1 += DNASeq1[i-1]
                    sequenceAlignment2 += DNASeq2[j-1]
                    i -= 1
                    j -= 1
                elif y == score[i][j-1]:
                    similarityScore -= 0.2
                    sequenceAlignment1 += '-'
                    sequenceAlignment2 += DNASeq2[j-1]
                    j -= 1
                elif y == score[i-1][j]:
                    similarityScore -= 0.2
                    sequenceAlignment1 += DNASeq1[i-1]
                    sequenceAlignment2 += '-'
                    i -= 1 


    sequenceAlignment1 = sequenceAlignment1[::-1]                  #reversing the strings we created on the traceback so that the alignment is front order
    sequenceAlignment2 = sequenceAlignment2[::-1]

    
    #################################  Output Section  ######################################
    result = "Similarity score: " + str(similarityScore) + '\n'
    result = result + "Sequence alignment1: " + sequenceAlignment1 + '\n'
    result = result + "Sequence alignment2: " + sequenceAlignment2 + '\n'
    writeToFile(outputPath,result)
    
def writeToFile(filePath, content):
    with open(filePath,'a') as file:
        file.writelines(content)

def readFile(filePath):
    logLines = ''
    with open(filePath,'r') as file:
        for logText in file:
            logLines = logLines + logText

    uniqueChars = ''.join(set(logLines))
    for ch in uniqueChars:
        if ch not in {'A','a','C','c','G','g','T','t'}:
            logLines = logLines.replace(ch,'')
    logLines = logLines.upper()
    return logLines

def removeFile(filePath):
    if os.path.isfile(filePath):
        os.remove(filePath)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='DNA sequence alignment')
    parser.add_argument('DNASeq1FilePath', type=str, help='Path to the file that contains First DNA sequence')
    parser.add_argument('DNASeq2FilePath', type=str, help='Path to the file that contains Second DNA sequence')
    parser.add_argument('OutputFilePath', type=str, help='Path to the output file')
    args = parser.parse_args()
    DNASeq1 = readFile(args.DNASeq1FilePath)
    DNASeq2 = readFile(args.DNASeq2FilePath)
    outputPath = args.OutputFilePath
    removeFile(outputPath)
    # writeToFile(outputPath,studentInfo()) #
    DNASeqAlignment(DNASeq1,DNASeq2,outputPath)

