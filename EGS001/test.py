mystring = "AANAA"
possiblematches = ["AAATA", "CCAAA", "AAAAC", "AAACA"]

for i, char in enumerate(mystring):
    if char == 'N':
        for match in possiblematches:
            corrected_string = mystring[:i] + match[i] + mystring[i+1:]
            if corrected_string in possiblematches:
                mystring = corrected_string
                break

print(mystring)
