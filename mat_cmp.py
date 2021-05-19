import sys

# Checking arguments
if len(sys.argv) != 3:
    print("No arguments provided, please provide two txt files.\n")
    sys.exit(1);

# Opening files
f1 = open(sys.argv[1], 'r')
f2 = open(sys.argv[2], 'r')

print("\nComparing files ", " @ refers to " + sys.argv[1], " # refers to " + sys.argv[2], sep='\n')

f1_line = f1.readline()
f2_line = f2.readline()

line_no = 1 # counter

with open(sys.argv[1]) as file1:
	with open(sys.argv[2]) as file2:
		same = set(file1).intersection(file2)

print("\nCommon Lines in Both Files")

for line in same:
	print(line, end='')

print('\n*****************************************************************************************')

print("\nDifferent Lines in Both Files")
while f1_line != '' or f2_line != '':
    # Removing whitespaces
    f1_line = f1_line.rstrip()
    f2_line = f2_line.rstrip()

	# Compares the lines from both file
    if f1_line != f2_line:
		
		# Outputs the line on file1 and use @ sign
        if f1_line == '':
            print("@", "Line-%d" % line_no, f1_line)
            print('---------------------------------------------------------------------')
        else:
            print("@-", "Line-%d" % line_no, f1_line)
            print('---------------------------------------------------------------------')
        
		# Outputs the line on file2 and use # sign
        if f2_line == '':
            print("#", "Line-%d" % line_no, f2_line)
            print('---------------------------------------------------------------------')
        else:
            print("#+", "Line-%d" % line_no, f2_line)
            print('---------------------------------------------------------------------')

        # Prints an empty line
        print()

	# Reads the next line from the file
    f1_line = f1.readline()
    f2_line = f2.readline()

    line_no += 1

f1.close()
f2.close()