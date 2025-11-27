package read;

import java.util.*;
import java.io.*;

/**
    Reads numerical constants from an input stream.

    Copyright 2012-2017 by Davide Bucci
    @author Davide Bucci
*/
public class CustomTokenizer
{
    private final Reader r;
    static private double[] powers10;

    /** Constructor.
        @param r the Reader to be used as the input stream.
    */
    public CustomTokenizer(final Reader r)
    {
        this.r = r;
    }

    /** Read a String up the end of the line.
    */
    public String nextLine() throws IOException
    {
        int ch;
        while ((ch=r.read())!='\n'&&ch!=-1) {
        }

        return "";
    }

    /** Read an int from the current stream (reader).
        This function "eats" all the spaces and characters which are preceding
        the double to be read. All the character of this list are ignored:
        ' ', '\r', '\t', ',', '(', ')'

        The following characters are interpreted as comments:
        '%', '#'
        and therefore the rest of the line will be ignored completely.

        @return the int constant read.
        @throws IOException when a problem occurs
    */
    public int nextInt() throws IOException
    {
        int i = handleIgnoredChars();
        if (i == -1)
            throw new EOFException();

        char c = (char) i;

        int result = (c - '0');
        while ((i = r.read()) >= 0) {
            c = (char) i;
            if (c == ' ' || c == '\n' || c == '\r' || (c - '0')>9 ||
                (c - '0')<0)
                break;
            result = result * 10 + (c - '0');
        }

        return result;
    }

    /** Ignore all characters which may be written before the double to be
        read from the stream. Handle comments as well.
        @return the last char read which is not an ignored character.
        @throws IOException when a problem occurs
    */
    private int handleIgnoredChars() throws IOException
    {
        int ch=r.read();
        while(ch==' '||ch=='\r'||ch=='\n'||ch=='\t'||ch==','||
                ch=='%'||ch=='#'||ch=='('||ch==')') {
            // '%' and '#' are comment characters: if we encounter them, we
            // skip all the line.
            if(ch=='%'||ch=='#'){
                while((ch=r.read())!='\n'&&ch!='\r'&&ch!=-1)
                    ;
                continue;
            } else {
                ch=r.read();
            }
        }
        return ch;
    }

    /** Read a double from the current stream (reader).
        This function "eats" all the spaces and characters which are preceding
        the double to be read. All the character of this list are ignored:
        ' ', '\r', '\t', ',', '(', ')'

        The following characters are interpreted as comments:
        '%', '#'
        and therefore the rest of the line will be ignored completely.

        @return the double constant read.
        @throws IOException when a problem occurs
    */
    public double nextDouble() throws IOException
    {
        final int NO_DECIMAL = -1;
        long number_value=0;
        double mantissa=0;
        boolean FlagScientific=false;
        boolean flagMinus=false;
        int iPointPos=NO_DECIMAL;
        int iCharCounter = 0;
        int ch;
        boolean signPossible = true;
        boolean  finish=false;


        ch=handleIgnoredChars();
        // We read a double-precision value here.
        while (!finish && (Character.isDigit(ch) ||ch==-1||
            ch=='.' || ch=='e' ||ch=='E'||ch=='-'||ch=='+')) {
            switch (ch) {
                case -1:
                    // Check wether a value has been read. signPossible is true
                    // only when the system expects to read a new value (which
                    // can be a mantissa or an exponent of a double precision
                    // constant. Therefore, if it is true, there is something
                    // missing and we raise an exception.
                    if(signPossible)
                        throw new EOFException();
                    finish=true;
                    break;
                case '-':
                    if (signPossible)
                        flagMinus = true;
                    // No break.
                case '+':
                    break;
                case '.':                   // Decimal separator.
                    iPointPos =iCharCounter+1;
                    break;
                case 'e':                   // Scientific-Notation 1E+14
                    // No break.
                case 'E':
                    FlagScientific=true;
                    // Store the value of the mantissa.
                    if (iPointPos!=NO_DECIMAL)     {
                        mantissa=number_value*pow10(iPointPos-iCharCounter);
                        iPointPos=NO_DECIMAL;
                    }else{
                        mantissa=number_value;
                    }
                    if(flagMinus)
                        mantissa=-mantissa;

                    flagMinus = false;
                    // Now number_value will contain the exponent.
                    number_value=0;
                    signPossible = true;
                    break;
                default:    // Process digits
                    signPossible = false;
                    number_value *= 10;
                    // If the value becomes negative, this means that we have
                    // an overflow here.
                    if(number_value<0)
                        throw new EOFException();
                    number_value += (int)(ch-'0');
            }
            ch=r.read();
            ++iCharCounter;
        }

        if (flagMinus)
            number_value= -number_value;  // Handle negative exponent

        double ff=0;
        if (!FlagScientific) {
            if (iPointPos!=NO_DECIMAL)
                ff = (double)number_value*pow10(iPointPos-iCharCounter);
            else
                ff = number_value;
        } else {
            if (iPointPos!=NO_DECIMAL)
                mantissa *= pow10(-iCharCounter+iPointPos);

            ff = mantissa*pow10((int)number_value);
        }

        if(iCharCounter==0)
            throw new EOFException();
        return ff;
    }

    /** Calculate (in a fast way) an integer power of 10.
        @param n the power to be calculated 10^n.
        @return the value of the power.
    */
    private static double pow10(int n)
    {
        double d=1.0;
        int k=n>0?n:-n;

        if(k<308) {
            if(powers10==null)
                powers10 = new double[308];

            if (powers10[k]==0) {
                for(int i=0;i<k;++i)
                    d*=10;
                powers10[k]=d;
            } else {
                d=powers10[k];
            }
        } else {
            d=Math.pow(10.0,n);
        }
        if(n<0)
            d = 1.0/d;

        return d;
    }
}

