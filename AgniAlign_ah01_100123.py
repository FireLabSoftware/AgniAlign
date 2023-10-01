#!/usr/bin/env python -i
## 
## ######################
## ->AgniAlign-- Really slow alignment for complex sequence relationships
##   *** Requires this script and the "VSG_ModuleFB.py" script in the same directory ***
##  -
## ->Version:
##  - This is version ah01 (Oct-01-2023).  This is an extreme pre-alpha version.  This is
##  - Provisional experimental software that seems to do something but might not behave as expected
##  - (which is not saying a lot, since it's not clear what is expected or what is "correct" in such alignments)
##
## ->Fuction: AgniAlign is used to provide a provisional picture of how a product from polymerase action
##  - relates to the initial template(s).  Of note, the program differs from conventional aligners in
##  - that the product is allowed to prime on itself or on the initial template in a manner that 
##  - allows complex products from a simple DNA or RNA input.  The program is intended for relatively short
##  - sequences (at most a few hundred bases) and is extremely slow and memory intensive for anything
##  - larger
##  - The algorithm here is a very highly modified Needleman Wunsch alignment with the "tweak" that
##  - both strands of the template are allowed as sources of alignment, along with (optionally) the
##  - product being allowed as a template for its own synthesis.  A major parameter choice is RefAlignOnly.
##  - Setting this to True (default) instructs the program to account for every base in the product.
##  - Untemplated bases are thus not allowed, but
##  - they will generally show up as isolated alignments to a specific position in the template
##  - Setting RefAlignOnly to False allows untemplated bases and may be more useful for many applications
##  - A template is not required, and setting Template='' will attempt to parse the product from its own sequence
##  -
##  ->BigPicture: AgniAlign is a parsing program for relatively short sequences, assigned either based on a template
##  - or based on assumptions that your molecule has served as its own template.
##  - The program is designed for relatively short sequences  (<300 bases) and will both given nonsensical results
##  - and take forever for longer sequences.
##  - In particular, AgniAlign is not a substitute for Blast, BWA, Bowtie, Star, Blat, or CandyCrush
##  - And likewise, the result is not a folding-- this is not a substitute for RNA fold
##  - The program will indicate some fundamentally complex repeated or internally complementary repetitive structures
##  - While the final parsing may not be unique, the may be instructive in understanding what might have occured during
##  - Formation of a molecule of interest
##  -  
## ->Syntax:
##  - python AgniAlign_<version> Template=<Sequence or FastA_FileName> Product=<Sequence_or_FastAFileName>
##  -    This attempts to parse a product from a defined template with no untemplated bases allowed
##  - python AgniAlign_<version> RefAlignOnly=False Template=<Sequence or FastA_FileName> Product=<Sequence_or_FastAFileName>
##  -    This attempts to parse a product from a defined template with untemplated bases allowed
##  - python AgniAlign_<version> RefAlignOnly=False Template='' Product=<Sequence_or_FastAFileName>
##  -    This attempts to parse a product from its own template
##  -
## ->Required Parameters
##  - Template=<Sequence or FastA_FileName> : Starting template sequence 
##  - Product=<Sequence or FastA_FileName>: List of product sequences (either in command line or from a FastA file)
##  -
## ->Optional Parameters
##  - Any of the valuoes in the "Parameters" section can be set from the command line with the syntax <parametername>=<Value>
##  -   There should be no spaces in these assignments (either in the value or around the eqauls sign
##  -
## ->Other Notes
##  - Requires Python 3.7+.  The pypy interpreter (www.pypy.org) may improve speed up to several fold.
##  - For Current Version of this code and the required VSG Modulce, see:  https://github.com/FireLabSoftware/AgniAlign
##  - Copywrite 2023 Andrew Fire and Stanford University, All Rights Reserved
##  - Inspiration/Advice: Emily Greenwald, Drew Galls, Orkan Ilbay, Dae Eun Jeong, Vanya Zheludev, David Lipman, Usman Enam, Karen Artiles, Nelson Hall, Matt McCoy, Janie Kim, Collete Benko, Tsachy Weissman, Nimit Jain, William Wang
## ###############
## End Help

## Begin List of Parameters
## Rewards and penalties used to choose the optimal alignment
## Note that these parameters will greatly affect the predicted alignments, particularly with no template
## Since first principles, there's no 'right' answer to either, this remains to user
## FWIW, the default parameters as of beta-release are listed for each 
Match1 = 10           ## (match=) Reward assigned if two bases match in a candidate alignment <default=10>
MisMatch1 = -20       ## (mismatch=) Cost assigned if two bases mismatch <default=-20>
GapStart1 = -40       ## (gapstart=) Cost assigned for starting a gap <default=-40> ## note that a gap here is a base in the template that is skipped over in the product-- not vice versa      
GapExtend1 = -2       ## (gapextend=) Cost assigned for extending a gap <default=-2>
Jump1 = -80           ## (jump=) Penalty for jumping between templates <default=-80>
SwitchStart1 = -40    ## (switchstart=) Base Penalty for switching between strands <default=-40>
SwitchExtend1 = -1    ## (switchextend=) Incremental (per base) Penalty for switching between strands <default=-1>
UnmatchedToMatched1 = -30 ## (unmatchedtomatched=) Penalty for switching from matched to unmatched <default=-30>
MatchedToUnmatched1 = -30 ## (matchedtounmatched=) Penalty for switching from unmatched to matched <default=-30>
HairpinTurn1 = -40    ## (hairpinturn=) Penalty for a zero-base-loop hairpin turn <default=-40>

## Algorithm Options
RefAlignOnly1 = True ## (refalignonly=) Setting this to False adds self-alignment capabilities
PrioritizeInitialMatch1 = 1 ## (priority=) For equivalent scores for two alignments etermines whether Squegee chooses the first match (1) or last (2)

## Display Paramaters
Font1 = "DejaVuSerifMono Bold" ## (font=) Font used for the output graphic
Scale1 = 24                ## (scale=) Scale for the drawing (and font size) pixels per base
OutputFile1 = 'default'    ## (outputfile=) Name of output file.  Can also be a svg, pdf, png, tiff, html, etc.  SVG is recommended
Rows1 = 'default' ## (rows=) Number of rows in the output (for multi-sequence figures)
TemplateRNAMode1 = 'Auto' ## (templaternamode=) 'Auto' sets an RNA mode using "U" bases instead of T if there are Us in the input template.  Setting to "True" forces this

## Default Sequences (defaults are included in code for Debug Purposes)
Product1 = ''## (product=) 'GGATTAAATTTCATATTGTTAATATTTATTAATGTATGTACAATATGAAATTTCATATTGTACATACATTAATAAATATTAACAATATGAAATTTCGG' ##(product=) Experimental sequence observed
Template1 = '' ## (template=) Template for the alignment (test for program development from Drew Galls was 'ACATACATTAATAAATATTAACAATATGAAATTTC')
ProductName1 = 'default' ## (productname=) A name for the product
TemplateName1 = 'default' ## (templatename=) A user name for the template
## for multiple products or templates, use commas to separate or a FastA input file

## Some colors that can be set
CircleColor1 = {('Input1',1):'blue', ('Input1',-1):'red', ('Output1',1):'blue',('Output1',-1):'red'} ## (circlecolor=) +1 are sense, -1 are antisense
TextColor1 = {True:'black', False:'gray50'} ## (textcolor=) Different text colors based on whether a given base is matched (true) or mismatched (false)

## End List of Parameters
import sys

if len(sys.argv)==2 and not('=' in sys.argv[1]) and not('h' in sys.argv[1].lower()) and all([x in 'agctuAGCTU' for x in sys.argv[1]]):
    Product1=sys.argv[1].strip().strip('"').strip("'").strip()
    sys.argv=sys.argv[:1]+sys.argv[2:]
from VSG_ModuleFB import *
from collections import Counter
vCommand()
if ProductName1=='default':
    ProductName1 = vnow+"_UserProduct"
if TemplateName1=='default':
    TemplateName1 = vnow+"_UserTemplate"
if Template1 == '':
    RefAlignOnly1 = False
if type(Product1)==str and os.path.isfile(Product1):
    ProductD1 = vFastAToDict(Product1)
else:
    if type(Product1)==str:
        Products1 = Product1.split(',')
    else:
        Products1 = Product1
    ProductNames1 = ProductName1.split(',')
    if len(ProductNames1)==1 and len(Products1)>1:
        ProductNames1 = [ProductName1+'_'+str(i+1) for i in range(len(Products1))]
    ProductD1 = dict(zip(ProductNames1,Products1))
if not ProductD1 or not(''.join(ProductD1.values())):
    print()
    print('******************************************************************')
    print('Warning: Agni not provide with a product to parse:')
    print("  Agni's Minimal syntax is")
    print("  python AgniAlign## Product=<Sequence_or_FastA_file>")
    print("  (with no spaces before or after equals in command line)")
    print('******************************************************************')
    print()
if type(Template1)==str and os.path.isfile(Template1):
    TemplateD1 = vFastAToDict(Template1)
else:
    if type(Template1)==str:
        Templates1 = Template1.split(',')
    else:
        Templates1 = Template1
    TemplateNames1 = TemplateName1.split(',')
    if len(TemplateNames1)==1 and len(Templates1)>1:
        TemplateNames1 = [TemplateName1+'_'+str(i+1) for i in range(len(Templates1))]
    TemplateD1 = dict(zip(TemplateNames1,Templates1))

if OutputFile1 == 'default':
    OutputFile1 = 'SqueegeeAlign'+'_'.join(ProductD1.keys())+'x'+'_'.join(TemplateD1.keys())+'_'+vnow+'.svg'.replace(' ','_').replace('@','a').replace(':','_')
myDefault1 = -99999999  ## An extremely low score used as a default

N1 = 'GATCRRRRKKKMMMMSSSSYYYYWWWWBBBBBVVVVVHHHHHHDDDDDDNNNNNNNUIJLFQEOOOOZFXZZTU-_=+P'
N2 = 'CTAGYCTUMACKTGUSCGIRAGIWATUVACGIBTCGUDATUGIHATCUINGATCIUACAEQFLCATUFZXTUZZ-_=+A'
AntiBaseS1 = set()
AllBases1 = set(list(N1))
AntiBaseD1 = {}
for n1,n2 in zip(N1,N2):
    AntiBaseS1.add((n1.upper(),n2.upper()))
    AntiBaseS1.add((n1.lower(),n2.upper()))
    AntiBaseS1.add((n1.lower(),n2.lower()))
    AntiBaseS1.add((n1.upper(),n2.lower()))
    AntiBaseS1.add((n2.upper(),n1.upper()))
    AntiBaseS1.add((n2.lower(),n1.upper()))
    AntiBaseS1.add((n2.lower(),n1.lower()))
    AntiBaseS1.add((n2.upper(),n1.lower()))
    if not(n1 in AntiBaseD1):
        AntiBaseD1[n1] = n2
    
def StripSeq1(s):
    sss = ''
    for b in s:
        if b.upper() in AllBases1:
            sss+=b
    return sss

def myArrow1(p1,p2,base,scale,zoner,width,color): ## This is just an algorithm that uses VSG to draw arrows
    if not(p1) or not(p2): return
    p1 = [z*scale for z in p1]
    p2 = [z*scale for z in p2]
    if p1==p2:
        return((p1,p1,p1))
    diffvect = (p2[0]-p1[0],p2[1]-p1[1])
    ldiffvect = (diffvect[0]**2+diffvect[1]**2)**0.5
    unitdiffvect = diffvect[0]/ldiffvect, diffvect[1]/ldiffvect
    orthounitvec = diffvect[1]/ldiffvect, -diffvect[0]/ldiffvect
    p1 = p1[0]+zoner*unitdiffvect[0], p1[1]+zoner*unitdiffvect[1]
    p2 = p2[0]-zoner*unitdiffvect[0], p2[1]-zoner*unitdiffvect[1]
    backarrowpoint = p2[0]-base*unitdiffvect[0], p2[1]-base*unitdiffvect[1]
    vline(x1=p1[0], y1=p1[1], x2=p2[0], y2=p2[1],stroke=color, strokewidth=width, priority=-1)
    vpolygon(points=((backarrowpoint[0]+base*orthounitvec[0]/2.0,backarrowpoint[1]+base*orthounitvec[1]/2.0),
           (backarrowpoint[0]-base*orthounitvec[0]/2.0,backarrowpoint[1]-base*orthounitvec[1]/2.0),
           (p2[0],p2[1])), fill=color, stroke=none, strokewidth=0, priority=-1)
def MatchScore1(o1,o2): ## returns the score for two position objects in terms of whether they match.  Not used currently by program
    if o1.b==o2.b:
        return Match1
    else:
        return MisMatch1
def ComplementScore1(o1,o2): ## returns the score for two position objects in terms of whether they are complementary
    if (o1.b,o2.b) in AntiBaseS1:
        return Match1
    elif o1.id==o2.id:
        return 0
    else:
        return MisMatch1
def ProximityScore1(o0,o1,o2,o3): ## returns a proximity score for two position object pairs in terms of whether they are adjacent and how far
    ## o0 and o1 are the previous and current positions in the product and these are (for the simple Aligner) consecutive in the product
    ## o3 and o2 are the previous and current templates.  They are by definition opposite strand unless o3==o0 and o1==o2 in which
    ## case these are unmapped to upstream sequence in product
    if o3.n != o2.n:
        return Jump1
    if o3.o != o2.o:
        return max(Jump1,SwitchStart1+SwitchExtend1*abs(o2.p-o3.p-o3.o))
    if o0.id==o3.id and o1.id==o2.id:
        return 0
    if o0.id==o3.id and o1.id!=o2.id:
        return UnmatchedToMatched1
    if o0.id!=o3.id and o1.id==o2.id:
        return MatchedToUnmatched1
    if o2.next==o3.id:
        return 0
    if o1.previous==o2.id:
        return HairpinTurn1
    else:
        return max(Jump1,GapStart1+GapExtend1*abs(o2.p-o3.p-o3.o))

class PositionObject1():
    def __init__(self,name,orientation,position,base,end):
        self.n = name ## An id name for sequence
        self.o = orientation ## 1 or -1
        self.oo = orientation ## will eventually note the orientation of the original template
        self.p = position ## zero based position in the sequence
        self.b = base ## base at that position
        self.id = (name,orientation,position)
        self.wc = (name,-orientation,position)
        self.t = None ## ID of template base
        self.s = 0 ## 0 based position of base in consecutive match stretch, numbered from 3' end
        self.v = 0 ## Score value for this position
        if position>0 and orientation==1:
            self.previous = (name,orientation,position-1)
        elif orientation==-1 and not(end):
            self.previous = (name,orientation,position+1)
        else:
            self.previous = None
        if not(end) and orientation==1:
            self.next = (name,orientation,position+1)
        elif orientation==-1 and position>0:
            self.next = (name,orientation,position-1)
        else:
            self.next = None
        self.x = None ## position on xy grid
        self.y = None ## position on xy grid
        self.e = end
        self.c = False ## True if complementary base to template; 0 if noncomplementary (mismatch)
        self.h = 'gray50' ## hue (color as an r,g,b tuple or a name
    def draw(self,poly):        
        vtext(text=self.b,
                  xc=self.x*Scale1,
                  yc=self.y*Scale1,
                  priority = 1,
                  font = Font1+' '+str(Scale1//2),
                  color = TextColor1[self.c])
        vtext(text=str(self.p+1),
                  xc=self.x*Scale1,
                  yc=self.y*Scale1-Scale1//4,
                  priority = 1,
                  font = Font1+' '+str(Scale1//5),
                  color = "gray50")
        if poly.lower().startswith('circ'):
            vcircle(xc=self.x*Scale1,
                yc=self.y*Scale1,
                r = Scale1//2.5,
                stroke = self.h,
                fill = white,
                strokewidth = 2,
                priority = 0)
        if poly.lower().startswith('squ') or poly.lower().startswith('rec'):
            vrectangle(xc=self.x*Scale1,
                yc=self.y*Scale1,
                r = Scale1//2.5,
                stroke = self.h,
                fill = white,
                strokewidth = 1,
                priority = 0)
        if self.previous:
            Previous=idD1[self.previous]
            myArrow1((Previous.x,Previous.y),(self.x,self.y),8,Scale1,Scale1/2.7,.8,black)
        if self.p==0 and self.o == 1:
            vtext(text=self.n,x1=self.x*Scale1-Scale1//4,yc=self.y*Scale1+Scale1,font=Font1+' '+str(Scale1),color="black") 

            

if (type(TemplateRNAMode1)==bool and TemplateRNAMode1==True) or (type(TemplateRNAMode1)==str and TemplateRNAMode1.lower()=='true'):
    UseU1 = True
else:
    UseU1 = False
if TemplateRNAMode1.lower()=='auto':
    if 'U' in ''.join(TemplateD1.values()):
        UseU1 = True
    
newy1 = 0
nextX0 = 0
if Rows1=='default':
    Rows1 = int((len(ProductD1)*3)**0.5)
for ProductNumber1,ProductName1 in enumerate(ProductD1):
    if ProductNumber1>0 and ProductNumber1%Rows1==0:
        nextX0 = int(VSG.xmax/Scale1)+3
        newy1 = 0
    nextX1 = nextX0
    idD0 = {}
    for n1 in TemplateD1:
        for p1,b1 in enumerate(StripSeq1(TemplateD1[n1])):
            end1 = False
            if p1+1==len(TemplateD1[n1]):
                end1 = True
            idD0[(n1,1,p1)] = PositionObject1(n1,1,p1,b1,end1)
            b11 = AntiBaseD1[b1]
            if UseU1 and b11.lower()=='T':
                b11 = b11.replace('T','U').replace('t','u')
            idD0[(n1,-1,p1)] = PositionObject1(n1,-1,p1,b11,end1)
    for poid1 in idD0:
        if idD0[poid1].p==0 and idD0[poid1].o == 1:
            nextX1 += 1
        idD0[poid1].x = nextX1
        idD0[poid1].y = 0
        idD0[poid1].h = CircleColor1[('Input1',idD0[poid1].o)]
        idD0[poid1].c = 1
        if idD0[poid1].o==-1:
            nextX1 += 1

    ThisProduct1 = ProductD1[ProductName1]
    idD1 = deepcopy(idD0) ## Individual positions in product or template. keys are ObjectIDs (name, orientation, position), values are PositionObjects
    ScoreD1 = {} ## keys are pairs of ids (reference, experimental), values are best score for that pair
    PreviousD1 = {} ## keys are pairs of ids (reference, experimental), values are best previous pair for that pair
    for p1,b1 in enumerate(StripSeq1(ThisProduct1)):
        end1 = False
        if p1+1==len(ThisProduct1):
            end1 = True
        idD1[(ProductName1,1,p1)] = PositionObject1(ProductName1,1,p1,b1,end1)
    endBestValue1 = myDefault1
    endBestCombination1 = None
    for poid1 in idD1:
        po1 = idD1[poid1]
        if not(po1.n==ProductName1):continue
        for poid2 in idD1:
            po2 = idD1[poid2]
            if po2.n==po1.n and (RefAlignOnly1 or po2.p>po1.p): continue ## enforces a requirement that a template base in the same molecule must be upstream
            if po1.p==0:
                ourBestValue1 = ComplementScore1(po1,po2) ## Start with a value of zero for each possible match
                ourBestPrevious1 = None
            else:
                ourBestValue1 = myDefault1
                ourBestPrevious1 = None
                for poid3 in idD1:
                    po3 = idD1[poid3]
                    if po3.n==po1.n and (RefAlignOnly1 or po3.p>=po1.p):
                        continue
                    po0 = idD1[po1.previous]
                    NewValue1 = ScoreD1[(poid3,po1.previous)]+ComplementScore1(po1,po2)+ProximityScore1(po0,po1,po2,po3)
                    if NewValue1>ourBestValue1 or (PrioritizeInitialMatch1==2 and NewValue1==ourBestValue1):
                        ourBestValue1 = NewValue1
                        ourBestPrevious1 = poid3,po1.previous
                if po1.e and (ourBestValue1>endBestValue1):
                    endBestValue1 = ourBestValue1
                    endBestCombination1 = (poid2,poid1)
            ScoreD1[(poid2,poid1)] = ourBestValue1
            PreviousD1[(poid2,poid1)] = ourBestPrevious1
    EChain1 = []
    nextCombo1 = endBestCombination1
    while nextCombo1:
        idD1[nextCombo1[1]].v = ScoreD1[nextCombo1]
        EChain1 = [nextCombo1,]+ EChain1
        nextCombo1 = PreviousD1[nextCombo1]
        
        
    lastpo2 = None
    xyD1 = {}
    for poid2,poid1 in EChain1:
        po1 = idD1[poid1]
        po2 = idD1[poid2]
        if (po1.b,po2.b) in AntiBaseS1:
            po1.c = True
        po1.t = po2.id
        if po2.oo==1:
            po1.oo = -1
            po1.h = CircleColor1['Output1',-1]
        else:
            po1.oo = 1
            po1.h = CircleColor1['Output1',1]
        BumpDown1 = False
        lastpo1 = None
        if po1.previous:
            lastpo1 = idD1[po1.previous]
        if po1.id==po2.id:
            if not(lastpo2) or lastpo1.id!=lastpo2.id:
                BumpDown1 = True
        else:
            if not(lastpo2):
                BumpDown1 = True
            elif po2.next and lastpo2.id!=po2.next:
                BumpDown1 = True
            elif lastpo2.oo!=po2.oo:
                BumpDown1 = True
            elif not((lastpo2.b,lastpo1.b) in AntiBaseS1):
                BumpDown1 = True
            elif not((po2.b,po1.b) in AntiBaseS1):
                BumpDown1 = True
        if (po1.x,newy1) in xyD1:
            BumpDown1 = True
        if po1.p==0:
            BumpDown1 = False
        if po1.id==po2.id:
            po1.x = nextX1
            nextX1 +=1
        else:
            po1.x = po2.x
        if (po1.x,newy1) in xyD1:
            BumpDown1 = True
        if BumpDown1:
            newy1 -= 1
        po1.y = newy1
        xyD1[(po1.x,newy1)] = po1
        po1.draw('circle')
        lastpo2 = po2
    newy1 -= 3
    for q1 in idD0:
        idD1[q1].y = newy1+(idD1[q1].o+1)//2
    for q1 in idD0:
        idD1[q1].draw('rectangle')
    newy1 -= 2
        
        

        
            
try:
    vdisplay(OutputFile1)
except:
    vdisplay(vnow+'_SqueegeeOutput.svg')

        


