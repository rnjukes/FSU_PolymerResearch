# Ryan Jukes
# 10/12/25
# Python Proficiency Quiz 1

#1 - Importing numpy and matplot into the file
import numpy as np
import matplotlib.pyplot as plt

#2 - Make list and textually define listed variables 
myFirstLast = [1,"dog",2.45,"apple",7.65]


for i in myFirstLast:
    if type(i) == str:
        print("this is a string")
    elif type(i) == int:
        print("this is an integer") 
    elif type(i) == float: 
        print("this is a float") 


#3 - Make the 5x5 array
A = np.zeros((5, 5))
print("Matrix A before edits:")
print(A)

#4 - Change the array to diagonal -2 and 1 pattern
i = 0
while i < 5:
    A[i, i] = -2
    if i > 0:
        A[i, i-1] = 1
    if i < 4:
        A[i, i+1] = 1
    i += 1  

print("Matrix A after edits:")
print(A)

print()


#5 - Make vector
b = np.ones((5,1))


#6 - Use linalg
x = np.linalg.solve(A,b)
print(x)


#7 - Make linspace array
print()
b = np.ones((5,1))
z = np.linspace(-np.pi,np.pi,100) 


#8 - Plot trig 
plt.plot(z,np.sin(z))
plt.plot(z,np.tanh(z))


#9 - Label trig plot
plt.xlabel("z")
plt.ylabel("f(z)")
plt.show() 
