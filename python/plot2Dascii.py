import sys

# Print a mode on the terminal.

def printmap(map):
    maxv=0
    ASCIIscale=[' ', '\033[37m.\033[0m','\033[33m+\033[0m',
        '\033[90m*\033[0m','\033[31mO\033[0m','\033[34mX\033[0m']

    for ll in map:
        for v in ll:
            if abs(v)>maxv:
                maxv=abs(v)

    def linepr(ll):
        print ("+", end=' ')
        for v in ll:
            sys.stdout.write("-")
        print ("+")
    
    linepr(ll);
    
    if maxv==0:
        maxv=0.1
    for ll in map:
        print ("|", end='')
        for v in ll:
            r=round(abs(v)/maxv*len(ASCIIscale))
            if r<len(ASCIIscale):
                sys.stdout.write( ASCIIscale[int(r)])
            else:
                sys.stdout.write("0")
        print ("|")

    linepr(ll);
    
    print ("Legend:")
    k=0
    for code in ASCIIscale:
        print (code + " = ", end=' ')
        t=float(k)/len(ASCIIscale)*maxv
        print("%.2f" % t, end=' ')
        k+=1;
    
    t=float(k)/len(ASCIIscale)*maxv
    print("0 = %.2f" % t, end='\n')
    
    linepr(ll)
