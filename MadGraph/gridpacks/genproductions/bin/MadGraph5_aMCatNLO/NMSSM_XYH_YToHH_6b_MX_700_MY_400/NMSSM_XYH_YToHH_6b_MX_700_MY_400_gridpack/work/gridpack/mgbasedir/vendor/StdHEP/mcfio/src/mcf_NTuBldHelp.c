/*******************************************************************************
*									       *
* mcf_NTuBldHelp.c --    						       *
*									       *
*	P. Lebrun, September 1995.					       *
*									       *
*******************************************************************************/
#include <stdio.h>
#include <Xm/Xm.h>
#include "help.h"
/*
** Menu lists
**
** Some of the menu list arrays have the number of items specified because
** the IBM C compiler complains: "Number of initializers cannot be greater
** than the number of aggregate members"
*/
                  
static helpMenuInfo

aboutDataModel = {"About this DDL", 'D', True, HELP_TEXT,
" This MOTIF application is a small GUI for a dedicated DDL (Data Definition \
Language) tool.  The Data Model underneath this language can be called the \
Generalized Ntuple Model, first implement in the context of HBOOK4 and PAW. \
We think that this data model is particularly well suited for Charm or B-Physics \
analysis.  Let us start with a simple example, where one has to study the \
the topology of a vertex, where the \
number of tracks in this vertex fluctuates from event to event. Thus, \
some information about this vertex is clearly of fixed size, such as the position of the \
vertex, it chi-squared and so forth, while some information has to be given \
on a track per track basis. So, in this model, we have some key variables \
to define first: the index pointing to a specific track, the current number \
of tracks for a particular instance of this d/s, called the Multiplicity, \
and the maximum possible value for this multiplicity across the entire \
event sample.  This last variable is defined given that, in F77, heap space \
for the d/s must be pre-allocated and can not be dynamically extended. \n\n\
	In our d/s, a given variable is therefore indexed, or, alternatively, \
is instanciated only once. In our example, we have only one value for the \
Vertex chi-square, while the track distance of closest approach to this \
vertex is a property of the track and is therefore indexed. Note that this \
type of data structure does not support abitrarily nesting of substructures: \
we simply stop at one level deep. Howver, within the context of MCFIO, \
one can build a hierarchical tree of multiple instances of such Generalized \
Ntuple Data structure, pretty much the same way one can define a tree for \
a set HistoScope Ntuple. \n\n\
	Although our Ntuple Data Model is quite simple, in VAX FORTRAN \
(or F90) \
and in C, it already has multiple possible implementations: one could treat \
the information about tracks (indexed variables) in parallell arrays: \n\n\
        INTEGER MAXNTR\n\
        PARAMETER (MAXNTR = 10)\n\
	COMMON/MYV/ NTR_IN_VERTEX, POSXYZ(3), CHI_V, \n\
     &     PARAMTR(10,MAXNTR), DIST_CLOSEST(MAXNTR) \n\n\
where PARAMTR, the track parameter of relevance in this study, and \
DIST_CLOSEST, the closest distance approach are multidimensioned array, with \
the slowest varying index being common to both variables, and is \
mapped to the track number index, running from 1 to MAXNTR. \
This type of COMMON has been \
extensively used in E687, and probably in other experiments as well. \n\n\
	However, in more modern languages, the following construct is a bit \
more rational:\n\n\
     structure /MYV_v_struct/ \n\
            Real            PARAMTR(10) \n\
            Character       DIST_CLOSEST \n\
      end structure\n\n\
     structure /MYV_struct/ \n\
            Real            NTR_IN_VERTEX \n\
            Real            POSXYZ(3) \n\
            Real            CHI_V     \n\
            record /MYV_v_struct/ track(nn2_max) \n\
     end structure\n\n\
Note that this DDL says nothing about internal organization of the data \
with respect to the event stream. One can group instances of Ntuple within \
events, keeping the information for one entry in the Ntuple together. \
In HBOOK4 terms, this is called the \"row-wise\" organization. Alternatively \
it became now fashionable to collect all entries for a particular variable \
together, allowing very fast histogramming of that variable (the \
 \"Column-wise\" organization).  Utilities could be written to convert \
one representation to the other, in the context of the MCFIO. \n\n",0,0},

aboutThisGUI = {"Why do I need a GUI?", 'P', False, HELP_TEXT,
" Why do I need to go through this moronic mouse exercise, instead of simply \
typing the COMMON block definition using my preferate editor? \n\n As stated \
above, you could type it two different ways, with no automatic way to change \
representations.  Moreover, it is likely that you'll find C and/or F90 \
customers for your data model, you would have to type the .inc and .h files \
separatly. Also, the straight F77 COMMON does not express unambiguously \
the array parallelism (e.g., which indices runs against the number of track \
in our example) while this information is saved in the file generated by this \
application.  Note also that variable casting is automatically inforced, \
character strings corresponding to variable names are checked for illegal \
characters, and so forth.. \n\n\
Finally, once your d/s is implemented using this \
tool, it can be declared to MCFIO in a straightforward and unambiguous \
fashion, and stored in an mcfio file so that data can be studied, histogrammed \
and displayed using the MCFIO data browser.  In other word, once a data model \
is formally expressed using this tool, the data written in this model context \
is fully self-described, and dedicated visualization tools\
 can be written. \n\n",0,0},
  

aboutThisPanel = {"About this Panel", 'P', False, HELP_TEXT,
" A Brief description on how to use this panel. \n\n\
This panel is split into two sections, the top part being dedicated to the \
global aspect of the generalized Ntuple to be created, the bottom part is \
used to define variables that are part of the Ntuple. In addition, a standard \
menu bar supports access to files, and allows to instanciated more  than \
one copy of this panel, and/or cloning data structures easily. \n\n\
Start by filling a title, or a name for the data structure you wish to \
create.  This is done by filling the text widget located at the top left. \
Also, the user is responsible for  maintaining the version token (top right \
text widget).  The multiplicity variable described above (NTR_IN_VERTEX \
in our example) should also be filled, as well the maximum value, which \
indirectly sets the size of the structure.  Underneath is a text field for a \
meaningful description of the Ntuple. You don't really have to fill that \
field, it's there for the quiche-eating software managers, who have nothin \
better to do than griping against hard working physicists like you. \n\n\
The second part of the panel allows you to fill the Ntuple data structure \
with variables. The variable list located at the middle of the panel \
can not be directly edited, however, selected a line allows you to  \
enter the variable properties, such as the name, type, dimensions \
and whether or not the variable is indexed by multiplicity. A comment line \
(maximum 71 characters) for each variable can optionally be filled. \
The set of four buttons on the top left allows you to copy this information \
 from/to an internal clipboard. This clipboard is shared among all panels, \
 alowing for easy merge of two distinct data structures. \n\n",0,0},
 
aboutVarOrdering = {"About Variable Ordering", 'P', False, HELP_TEXT,
" Feel free to enter the variable in any order you whish. You can also \
modify this list in random order, or leave blank entry in this list. However \
be aware that, upon saving this information for future references, \
this application will reordered \
this list following this simple algorythm: first we place all the fixed size \
variable (e.g., non indexed), than the indexed variables. Within such \
subgroups, the order follows the list displayed in the pull-down menu option\
 from which  the variable type has been selected: first the BYTE, last the \
DOUBLE COMPLEX. \n\n\
Such is re-ordering allows for optimisation of the XDR filitering unit: \
By bunching variables of the same type, on many architecture, space is \
allocated without padding, (especially for the commonly used REAL, DOUBLE \
PRECISION), allowing efficient use of the XDR array function.  In addition, \
this might a fairly natural organization of the data. \n\n\
After this re-ordering is done, this application also check for \
alignements within what could be instanciated as F77 COMMON. To avoid \
confusion when such a COMMON is directly accessed, on many systems, it is \
necessary to align variables on a word boundary. The algorithm is actually \
simple to express mathematically. Let START be the starting address of the \
COMMON. Then one must have for each variable in the structure: \n\n\
  Modulus((Variable_address - start),(size of variable))=0 \n\n\
For instance, the following sequence is misaligned: \n\n\
               CHARACTER*5 NAME\n\
               DOUBLE PRECISION VALUE\n\n\
Obviously, it is safer to allocate 8 bytes for the character array, \
and not 5.  This application will let you place padding variables where \
deemed necessary, you can use them, or, conversely you can fix one of the \
array sizes.  Whatever is most convenient for your application, \
provided the above equation is respected.  Note that if VAX FORTRAN \
d/s are used, the variable within a substructure must also be aligned \
with respect to the starting point of the substructure. \n" , 0, 0}, 
 
inputFiles = {"dBin Input Files", 'I', False, HELP_TEXT,
"There is only one type of input file for this application, that is, files \
that have been created by this application, written using the dbin utility, \
written by T. Wenaus...\n",0,0},
                   
outputFiles = {"Output Files", 'O', False, HELP_TEXT,
" Brief description of the output files generated by this application. \n\n\
1. DDl Files \n\
The DDL for your particular data model can be saved as a dBin file. \
Such a dbin file is a small ASCII files, describing the NTuple \
and related variables.  It is probably a bad idea to temper with it, as it may \
become undreadable by this application. The pull-Down menu \"Save\" or \"Save as...\" \
triggers the creation or updating of this output file.  Note that, in order \
to create other types of files, you must first saved the DDL file, as the \
other files directly derived from this DDL file. \n\n\
2. The F77 or VAX FORTRAN Include file\n\
Directly Derived from the DDl content, including the name of this include \
file: if the DDL file name is xyz.db, this file will be xyz.inc. Depending \
on the chosen layout (see preference menu, second option), either a F77 \
(albeit with variable names that can be 31 char. long, and may contain \
byte arrays, and int*2) COMMON block statement is issued, or alternatively, \
a set of two nested structures VAX FORTRAN structures, and a COMMON \
containing a RECORD holding the top level structure.\n\n\
3. The ANSI C include file.\n\
Directly Derived from the DDl content, including the name of this include \
file: if the DDL file name is xyz.db, this file will be xyz.h. Depending on \
the chosen layout (see preference menu, second option), a single struture is \
defined, holding parallel arrays, or, alternatively, a struture holding the \
indexed variables, and a top level structure holding the fixed size \
quantities and a fixed number of child structures. See example above. \n",0,0},

preferences = {"Preferences", 'P', False, HELP_TEXT,
" Brief descriptions of various modes or preferences. \n\n\
1.Set Lang. env...\n\
Syntax for variable types, array ordering and index range is obviously \
different in ANSI C (or C++) , F77 and F90.  Currently only two modes are \
supported: F77 and C. Default is F77. \n\n\
2. Set Data Struct. Organisation \n\
A variable within the Ntuple comes as indexed, or not indexed. One can choose \
to group indexed variables in a sub-structure, or conversely, implement them \
as parallel arrays. See example above. Default is parallel arrays. \n",0,0},

terms = {"Terms", 'T', False, HELP_TEXT,
"Brief descrption of technical term used throughout this application:\n\n\
DDL: Data Definition Language: The language that help us implement a particular\
Data Descriptor, within the context of the Generalized Ntuple Data model.\n\n\
Version: On the panel, upper right hand corner, a 7-character long field is \
reserved for a user-defined, user-maintained token that states the particular \
version of his Data Descriptor. This token will be checked against the value \
found at run-time within MCFIO.\n",0,0},

relNts = {"Release Notes", 'R', False, HELP_TEXT,
"None for now.. \n",0,0},

vR = {"Version", 'V', False, HELP_TEXT,
"nTuple Data Definition GUI, Version 1.0\n\
Copyright (c) 1995, 1996 Universities Research Association, Inc.\n\
All rights reserved.\n\
\n\
This is an ancillary product for the mcfio system, the i/o system for McFast. \
written by the Simulation Group, reviewed by the members of the Nirvana\
 project. Some MOTIF utilities have been \
borrowed from the Fermilab Nirvana project.\n", 0,0};

helpMenuInfo

*NTuBldMenuHelp[]={&aboutDataModel, &aboutThisGUI, &aboutThisPanel,
                  &aboutVarOrdering, &inputFiles, &outputFiles,
                  &preferences,  &terms, &relNts, &vR, NULL};
