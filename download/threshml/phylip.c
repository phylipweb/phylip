/* version 3.7. (c) Copyright 1993-2005 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, Andrew Keeffe,
   Dan Fineman, Patrick Colacurcio, and Mike Palczewski.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#ifdef OSX_CARBON
#include <Carbon/Carbon.h>
#endif /* OSX_CARBON */

#include <stdio.h>
#include <signal.h>
#include "phylip.h"

#ifdef WIN32
#include <windows.h>
/* for console code (clear screen, text color settings) */
CONSOLE_SCREEN_BUFFER_INFO      savecsbi;
boolean    savecsbi_valid = false;
HANDLE  hConsoleOutput;

void phyClearScreen();
void phySaveConsoleAttributes();
void phySetConsoleAttributes();
void phyRestoreConsoleAttributes();
void phyFillScreenColor();
#endif /* WIN32 */

#include "phylip.h"
#include "Slist.h"

#ifndef OLDC
static void crash_handler(int signum);
boolean unrooted_tree_locrearrange_recurs(tree* t, node *p, node*pp, 
    double* bestyet, boolean thorough, tree* priortree, tree* bestree);
static void rooted_repreorder(tree* t, node *p, boolean *success);
static void rooted_tryrearr(tree*t, node *p, boolean *success);
static void _fgetline_finalize(void);
#endif /* OLDC */

#if defined(OSX_CARBON) && defined(__MWERKS__)
boolean fixedpath = false;
#endif /* OSX_CARBON && __MWERKS__ */

/* Global file objects */
/* TODO Use locals and/or move to individual programs */
FILE *infile, *outfile, *intree, *intree2, *outtree;
FILE *weightfile, *catfile, *ancfile, *mixfile, *factfile;

long spp;                       /* number of species */
long words, bits;
boolean ibmpc, ansi, tranvsp;
naym *nayme;                    /* names of species */
initdata functions;             /* pointers to specialized node_new
                                   and tree_new */

/* Default vtable for generic nodes. */
struct node_vtable node_vtable = {
  generic_node_init,
  generic_free_node,
  generic_node_copy
};


void no_op(void)
{ /* Do nothing. Used as a dummy pointer to a function that */
} /* doesn't need to do anything (e.g. smooth for parsimony)*/

/********* Tree and node functions ***********/

void even_sibs(tree* t, node* src, node* dst)
{ /* Add or delete nodes in dst fork until it has the same number as src
   * fork. New nodes are added or removed at dst->next.
   * src and dst must be fork nodes in t.
   *
   * Nothing is done to free or initialize adjacent tips or forks.
   */
  
  long src_sibs, dst_sibs;
  node* n;
  
  src_sibs = count_sibs(src);
  dst_sibs = count_sibs(dst);

  while ( src_sibs > dst_sibs) {
    n = t->get_forknode(t, src->index);
    n->next = dst->next;
    dst->next = n;
    dst_sibs++; 
 }
  
  while ( dst_sibs > src_sibs) {
    n = dst->next;
    dst->next = dst->next->next;
    t->release_forknode(t, n);
    dst_sibs--;
  }
}


node* where_in_dest (tree* src, tree* dst, node* nsrc )
{ /* Return the node in dst that corresponds to node nsrc in src
   * i.e. has the same index and is the same number of next links
   * away from src's nodep entry.
   *
   * Corresponding forks in src and dst should have equal number
   * of nodes beforehand.
   */
  
  node *ret = NULL, *p;
 
  if ( nsrc ) {
    p = src->nodep[nsrc->index - 1];
    ret = dst->nodep[nsrc->index - 1];
    while ( p != nsrc ) {
      p = p->next;
      ret = ret->next;
    }
  } 
  return ret;
}


void generic_tree_copy(tree* src, tree* dst)
{ /* copies tree src to tree dst*/
  long i, j, num_sibs;
  node *p, *q;
  Slist_node listnode;

  /* ensure that all forks in dst have same number of sibs as in src */
  for ( i = spp ; i < src->nonodes ; i++ ) 
    even_sibs(dst, src->nodep[i], dst->nodep[i]);

  /* copy tip nodes and link to proper dst forks */
  for (i = 0; i < spp; i++) {
    src->nodep[i]->copy(src->nodep[i], dst->nodep[i]);
    dst->nodep[i]->back = where_in_dest(src, dst, src->nodep[i]->back);
  }

  /* copy fork nodes and back links */
  for (i = spp; i < src->nonodes; i++) {
    p = src->nodep[i];
    q = dst->nodep[i];
    
    num_sibs = count_sibs(p);
    
    for (j = 0; j <= num_sibs; j++) { 
      p->copy(p, q);
      q->back = where_in_dest(src, dst, p->back);
      p = p->next;
      q = q->next;
    }
  }

  /* copy score and root */
  dst->score = src->score;
  dst->root = where_in_dest(src, dst, src->root);

  /* clear dst's list of free forks */
  /* Forks are still tracked in nodep */
  while ( !Slist_isempty(dst->free_forks) )
    Slist_pop(dst->free_forks);

  /* Iterate through src->free_forks, adding corresponding fork node to dst->free_forks. */
  for ( listnode = src->free_forks->first;
        listnode != NULL;
        listnode = listnode->next
      ) {
    /* Find the corresponding (i.e. linked by nodep) node in dst */
    p = where_in_dest(src, dst, listnode->data);
    Slist_append(dst->free_forks, p);
  }
}


void generic_node_copy(node* src, node* dst) 
{
  /* Copy node data from src to dst.
   *
   * FIXME how do we want this to work?  we probably don't want to copy
   * next and back, but copy everything else 
   * some already needed things are here
   */
  dst->v = src->v;
  dst->xcoord = src->xcoord;
  dst->ycoord = src->ycoord;
  dst->ymin = src->ymin;
  dst->ymax = src->ymax;
  dst->iter = src->iter;
  dst->haslength = src->haslength;
  dst->initialized = src->initialized;
  dst->deltav = src->deltav;

}

void generic_fork_print(node * n)
{
    boolean firstTime = true;
    boolean nulledOut = false;
    node * p = n;
    while((firstTime || (p != n)) && !nulledOut)
    {
        if(p == NULL)
        {
            nulledOut = true;
        }
        else
        {
            printf("  ");
            p->node_print_f(p);
            p = p->next;
            printf("\n");
        }
        firstTime = false;
    }
}

void generic_node_print(node *n)
{
    printf("%10p : %10p",n,n->back);
    if(n->back != NULL)
    {
        printf(" %3ld",n->back->index-1);
    }
    else
    {
        printf("    ");
    }
    printf(" p->v : %lf",n->v);
    printf(" p->iter : %d",n->iter);
    printf(" init : %d",n->initialized);
}


void generic_free_node(node **n) 
{ /* Release a node's memory */
  free(*n);
  *n = NULL;
}


/* generic_node_init
 *
 * Assign default node data. tip is set false when type is FORK_NODE (0)
 * otherwise true. Index is assigned as given.
 */
void generic_node_init(node* n, node_type type, long index)
{
  if ( type == TIP_NODE )
    n->tip = true;
  else if ( type == FORK_NODE )
    n->tip = false;
  else
    /* Since we used to simply pass type=true for tips, any other value
     * will be a tip... for now. */
    n->tip = true;

  n->index = index;
  n->v = initialv;
  n->iter = true;
  n->initialized = false;

  /* Initialize virtual functions */
  n->init = generic_node_init;
  n->free = generic_free_node;
  n->copy = generic_node_copy;
  n->reinit = generic_node_reinit;
  n->fork_print_f = generic_fork_print;
  n->node_print_f = generic_node_print;
}


void generic_node_reinit(node * n)
{
  n->back = NULL;
  n->v = initialv;
  n->iter = true;
  n->initialized = false;
  /*
     don't change n->index !!!
     code relies upon having a unique index for each node
     whether or not it's currently in the tree
  */
}


node* generic_new_node(node_type type, long index) 
{ /* Allocate, initialize, and return a new node, setting tip and index. */
  node* n = Malloc(sizeof(node));
  generic_node_init(n, type, index); 
  return n;
}

void gnu(node **grbg, node **p)
{ /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */

  if (*grbg != NULL) {
    *p = *grbg;
    *grbg = (*grbg)->next;
  } else
    *p = functions.new_node(false, 0);

  (*p)->back       = NULL;
  (*p)->next       = NULL;
  (*p)->init(*p, false, 0);
}  /* gnu */


void chuck(node **grbg, node *p)
{ /* collect garbage on p -- put it on front of garbage list */
  p->back = NULL;
  p->next = *grbg;
  *grbg = p;
}  /* chuck */


void chucktreenode(node **grbg, node *p)
{ /* collect garbage on p -- put it on front of garbage list */

  p->back = NULL;
  p->next = *grbg;
  *grbg = p;
}  /* chucktreenode */


void setupnode(node *p, long i)
{ /* initialization of node pointers, variables */

  p->next = NULL;
  p->back = NULL;
  p->index = i;
  p->tip = false;
}  /* setupnode */


long count_sibs (node *p)
{ /* Count the number of nodes in a ring, return the total number of */
  /* nodes excluding the one passed into the function (siblings)     */
  node *q;
  long return_int = 0;

  if (p->tip) {
   printf ("Error: the function count_sibs called on a tip.  This is a bug.\n");
    exxit (-1);
  }

  q = p->next;
  while (q != p) {
    if (q == NULL) {
      printf ("Error: a loop of nodes was not closed.\n");
      exxit (-1);
    } else {
      return_int++;
      q = q->next;
    }
  }

  return return_int;
}  /* count_sibs */


void verify_nuview(node *p)
{ /* DEBUG function. Traverses entire tree and prints error message
   * if any view towards p has not been initialized. */

  /* TODO: implement */
}

void invalidate_nuview(node *p)
{ /* Invalidate all views looking toward p. Must be called on a node
   * after changing its tyme or branch lengths before evaluating at any other
   * node. */
  
  invalidate_traverse(p);
  invalidate_traverse(p->back);
}

void invalidate_traverse(node *p)
{ /* Invalidates p's view and all views looking toward p from p->back
   * on out. */
  node *q;
  
  if (p == NULL)
    return;
  if (p->tip)
    return;

  p->initialized = false;

  q = p->back;
  if ( q == NULL ) return;
  if ( q->tip ) return;

  /* Call ourselves on p->back's sibs */
  for ( q = q->next ; q != p->back ; q = q->next) {
    invalidate_traverse(q);
  }
} /* invalidate_traverse */


void inittrav_all(tree *t)
{
  /* Set initialized false on all nodes reachable from p, so
   * that view are regenerated regardless. For debugging nuview
   * problems. */

  node *p;
  long index;

  /* For each fork node (spp..nonodes-1) */
  for ( index = spp; index < t->nonodes; index++ ) {
    p = t->nodep[index];
    
    /* Set initialized false on all nodes in fork */
    p->initialized = false;
    for ( p = p->next; p != t->nodep[index]; p = p->next ) {
      p->initialized = false;
    }
  }
} /* inittrav_all */


void inittrav (node *p)
{ /* traverse to set pointers uninitialized on inserting */
  node *sib_ptr;

  /*
  printf ("EWFIX.INIT %p\n",p);
  */
  
  if (p == NULL)
    return;
  if (p->tip)
    return;
  for ( sib_ptr  = p->next ; sib_ptr != p ; sib_ptr = sib_ptr->next) {
    /*
    printf("EWFIX -- uninit %p\n",sib_ptr);
    */
    sib_ptr->initialized = false;
    inittrav(sib_ptr->back);
  }
} /* inittrav */


/********* Error handling ***********/

void EOF_error(void)
{ /* Print a message and exit when EOF is reached prematurely. */
  puts("\n\nERROR: Unexpected end-of-file.\n");
  exxit(-1);
}


static void crash_handler(int sig_num)
{ /* If (when?) we crash, print out something useful */
  boolean segorbus;

  printf("ERROR:  ");
  switch(sig_num) {
#ifdef SIGSEGV
    case SIGSEGV:
      puts("This program has caused a Segmentation fault.");
      break;
#endif /* SIGSEGV */
#ifdef SIGFPE
    case SIGFPE:
      puts("This program has caused a Floating Point Exception");
      break;
#endif  /* SIGFPE */
#ifdef SIGILL
    case SIGILL:
      puts("This program has attempted an illegal instruction");
      break;
#endif  /* SIGILL */
#ifdef SIGPIPE 
    case SIGPIPE:
      puts("This program tried to write to a broken pipe");
      break;
#endif  /* SIGPIPE */
#ifdef SIGBUS
    case SIGBUS:
      puts("This program had a bus error");
      break;
#endif /* SIGBUS */
  }   
  segorbus = false;
#ifdef SIGBUS
  segorbus = segorbus || SIGBUS;
#endif /* SIGBUS */
#ifdef SIGSEGV
  segorbus = segorbus || SIGSEGV;
#endif /* SIGSEGV */
  if (segorbus) {
    puts(
    "       This may have been caused by an incorrectly formatted input file");
    puts(
        "       or input tree file.  You should check those files carefully.");
    puts("       If this seems to be a bug, please mail joe@gs.washington.edu");
  }
  else {
    puts("       Most likely, you have encountered a bug in the program.");
    puts("       Since this seems to be a bug, please mail joe@gs.washington.edu");
  }
  puts("       with the name of the program, your computer system type,");
  puts("       a full description of the problem, and with the input data file.");

#ifdef WIN32
  puts ("Press Enter or Return to close program.");
  puts("  You may have to press Enter or Return twice.");
  fgetline(stdin);
#endif

#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  abort();
}

/********** Initialization *************/

void init(int argc, char** argv, initdata* ini) 
{ /* initialization routine for all programs 
   * anything done at the beginning for every program should be done here */ 
 
  /* set up signal handler for segfault, floating point exception, illegal
   * instruction, bad pipe, bus error.  There are more signals that can cause a
   * crash, but these are the most common even these aren't found on all
   * machines.  */

#ifdef SIGSEGV
  signal(SIGSEGV, crash_handler);
#endif /* SIGSEGV */
#ifdef SIGFPE
  signal(SIGFPE, crash_handler);
#endif /* SIGFPE */
#ifdef SIGILL
  signal(SIGILL, crash_handler);
#endif /* SIGILL */
#ifdef SIGPIPE
  signal(SIGPIPE, crash_handler);
#endif /* SIGPIPE */
#ifdef SIGBUS
  signal(SIGBUS, crash_handler);
#endif /* SIGBUS */

  /* Set default terminal characteristics */
  ibmpc = IBMCRT;
  ansi = ANSICRT;

  /* Clear the screen */
  cleerhome();

  /* Perform DOS console configuration */
  phySetConsoleAttributes();
  phyClearScreen();

  /* initialize 'functions' as given, or provide defaults */
  if ( ini == NULL ) {
    functions.new_node = generic_new_node; 
    functions.tree_new = generic_tree_new;
  } else {
    if ( ini->new_node != NULL )
      functions.new_node = ini->new_node; 
    else 
      functions.new_node = generic_new_node; 
    if ( ini->tree_new != NULL )
      functions.tree_new = ini->tree_new;
    else 
      functions.tree_new = generic_tree_new;
  }
}

/************* File reading *************/

void scan_eoln(FILE *f) 
{ /* Eat everything up to EOF or newline, including newline */
  char ch;

  while (!eoff(f) && !eoln(f)) 
    gettc(f);
  if (!eoff(f)) 
    ch = gettc(f);
}


boolean eoff(FILE *f)
{ /* Return true iff next getc() is EOF */
    int ch;

    if (feof(f)) 
      return true;
    ch = getc(f);
    if (ch == EOF) {
      ungetc(ch, f);
      return true;
    }
    ungetc(ch, f);
    return false;
}  /*eoff*/


boolean eoln(FILE *f)
{ /* Return true iff next getc() is EOL or EOF */
    register int ch;

    ch = getc(f);
    if (ch == EOF)
      return true;
    ungetc(ch, f);
    return ((ch == '\n') || (ch == '\r'));
}  /*eoln*/


boolean filexists(const char *filename)
{ /* Return true iff file already exists */
  FILE *fp;
  fp = fopen(filename,"r");
  if (fp) {
    fclose(fp);
    return 1;
  } else
    return 0;
}  /*filexists*/


void openfile(
/* Attempt to open a file.
 * 
 * If file cannot be opened or already exists, the user is asked what to do.
 */
    FILE **fp,                  /* where to return FILE* */
    const char *filename,       /* file name to open */
    const char *filedesc,       /* description of file ("Input tree file") */
    const char *mode,           /* access mode to attempt */
    const char *application,    /* name of current program */
    char *perm                  /* buffer to return actual access permissions
                                   granted (may be NULL) */
)
{
  FILE *of;
  char file[FNMLNGTH];
  char filemode[3];
  char ch;
  const char *progname_without_path;
  long loopcount, loopcount2;
#if defined(OSX_CARBON) && defined(__MWERKS__)
  ProcessSerialNumber myProcess;
  FSRef bundleLocation;
  unsigned char bundlePath[FNMLNGTH];

  if(!fixedpath){
    /* change path to the bundle location instead of root directory */
    GetCurrentProcess(&myProcess);
    GetProcessBundleLocation(&myProcess, &bundleLocation);
    FSRefMakePath(&bundleLocation, bundlePath, FNMLNGTH);
    chdir((const char*)bundlePath);
    chdir(".."); /* get out of the .app directory */
    
    fixedpath = true;
  }
#endif

  progname_without_path = get_command_name(application);

  strcpy(file, filename);
  strcpy(filemode, mode);
  loopcount = 0;
  while (1) {
    /* Ask before clobbering existing file */
    if (filemode[0] == 'w' && filexists(file)) {
      printf("\n%s: the file \"%s\" that you wanted to\n",
          progname_without_path, file);
      printf("     use as %s already exists.\n", filedesc);
      printf("     Do you want to Replace it, Append to it,\n");
      printf("     write to a new File, or Quit?\n");
      loopcount2 = 0;
      for (;;) {
        printf("     (please type R, A, F, or Q) \n");
        phyFillScreenColor();
        fflush(stdout);
        if ( (ch = menu_getchar()) != '\0' ) {
          if (strchr("RAFQ", ch) != NULL )
            break;
        }
        countup(&loopcount2, 10);
      }
      switch (ch) {
        case 'R':
          /* do nothing special */
          break;
        case 'A':
          strcpy(filemode,"a");
          continue;
        case 'F':
          file[0] = '\0';
          loopcount2 = 0;
          while (file[0] =='\0') {
            printf("Please enter a new file name> ");
            fflush(stdout);
            getstryng(file);
            countup(&loopcount2, 10);
          }
          strcpy(filemode,"w");
          continue;
        case 'Q':
          exxit(-1);
          break;
        default:        /* Shouldn't happen */
          assert(0); 
      }
    }

    /* Open the file */
    of = fopen(file,filemode);
    if (of)
      break;
    else {
      /* Can't open. use mode flags to guess why not. */
      switch (filemode[0]){

      case 'r':
        printf("%s: can't find %s \"%s\"\n", progname_without_path,
            filedesc, file);
        file[0] = '\0';
        loopcount2 = 0;
        while ( file[0] =='\0' ) {
          printf("Please enter a new file name> ");
          fflush(stdout);
          getstryng(file);
          countup(&loopcount2, 10);
        }
        break;

      case 'w':
      case 'a':
        printf("%s: can't write %s \"%s\"\n", progname_without_path,
            filedesc, file);
        file[0] = '\0';
        loopcount2 = 0;
        while (file[0] == '\0') {
          printf("Please enter a new file name> ");
          fflush(stdout);
          getstryng(file);
          countup(&loopcount2, 10);
        }
        continue;
      default:
        printf("Internal error in openfile(). Unknown mode \"%s\"\n", filemode);
        exxit(-1);
      }
    }
    countup(&loopcount, 20);
  }
  *fp = of;
  if (perm != NULL)
    strcpy(perm, file);
} /* openfile */


const char* get_command_name (const char *vektor)
{ /* returns the name of the program from vektor without the whole path */
  char *last_slash;

  /* Point to the last slash... */
  last_slash = strrchr (vektor, DELIMITER);

  if (last_slash)
    /* If there was a last slash, return the character after it */
    return last_slash + 1;
  else
    /* If not, return the vector */
    return vektor;

}  /*get_command_name*/

/****************** User input ************/

/* Only fgetline() and fgetline_finalize() may use this. */
static char *_fgetline_buffer = NULL;

/* Only fgetline() may use this. */
static void _fgetline_finalize()
{
  /* Free dynamic memory used by fgetline */

  if (_fgetline_buffer != NULL) {
    free(_fgetline_buffer);
    _fgetline_buffer = NULL;
  }
}


char *fgetline(FILE *fp)
{
  /* Read a complete line of input, strip newline, and return a pointer to the
   * buffer. If EOF is encountered, program is aborted. Return value should be
   * considered read-only and not valid across calls to this function. */

  static size_t size = 0x100;
  char *b = NULL;
  char *lastch;
  unsigned long len;
  int i;

  assert(fp);

  /* Set function to free buffer on exit */
  if ( _fgetline_buffer == NULL ) {
    _fgetline_buffer = malloc(size);
    i = atexit(_fgetline_finalize);
    assert(i == 0);
  }

  b = _fgetline_buffer;
  for (;;) {
    if ( fgets(b, size - (b - _fgetline_buffer), fp) == NULL )
      EOF_error();
    len = strlen(_fgetline_buffer);
    lastch = _fgetline_buffer + len - 1;
    /* Check for newline, chomp, return */
    if ( *lastch == '\n' ) {
      (*lastch) = '\0';
      return _fgetline_buffer;
    }
    else {
      /* Double the size of buffer and continue */
      size *= 2;
      _fgetline_buffer = realloc(_fgetline_buffer, size);
      b = _fgetline_buffer + len;
    }
  }
}


char menu_getchar()
{
  /* Read a line from stdin and returns the uppercase version of the first
   * non-whitespace char, or '\0' if no such char exists. Aborts if EOF is
   * reached.  */
  
  char *line;
  char ch;
  int result;

  line = fgetline(stdin);       /* abort on EOF */
  result = sscanf(line, " %c", &ch);
  if ( result == 1 )
    return toupper(ch);

  return '\0';
}


void getstryng(char *fname)
{ /* read in a file name from stdin and take off newline if any */
  char *end;

  assert(fname);

  fname = fgets(fname, 100, stdin);
  if ( fname == NULL )
    EOF_error();

  if ( (end = strpbrk(fname, "\n\r")) != NULL)
    *end = '\0';
    
} /* getstryng */


void countup(long *loopcount, long maxcount)
{ /* count how many times this loop has tried to read data, bail out
     if exceeds maxcount */

  assert(loopcount);
  assert((*loopcount) < maxcount);

  (*loopcount)++;
  if ((*loopcount) >= maxcount) {
    printf("\nERROR: Made %ld attempts to read input in loop. Aborting run.\n",
            *loopcount);
    exxit(-1);
  }
} /* countup */


void cleerhome()
{ /* home cursor and clear screen, if possible */
#ifdef WIN32
  if(ibmpc || ansi){
    phyClearScreen();
  } else {
    printf("\n\n");
  }
#else
  printf("%s", ((ibmpc || ansi) ? ("\033[2J\033[H") : "\n\n"));
#endif
} /* cleerhome */


long readlong(const char *prompt)
{ 
  /* read a long */

  long res, loopcount;
  char *string;

  loopcount = 0;
  for (;;) {
    printf("%s", prompt);
    fflush(stdout);
    string = fgetline(stdin);
    if (sscanf(string,"%ld",&res) == 1)
      break;
    countup(&loopcount, 10);
  }

  return res;
}  /* readlong */


void uppercase(Char *ch)
{ /* convert ch to upper case */
  *ch = (islower (*ch) ? toupper(*ch) : (*ch));
}  /* uppercase */


/**************  Random number generation *********/

double randum(longer seed)
{ /* random number generator -- slow but machine independent
  This is a multiplicative congruential 32-bit generator
  x(t+1) = 1664525 * x(t) mod 2^32, one that passes the
  Coveyou-Macpherson and Lehmer tests, see Knuth ACP vol. 2
  We here implement it representing each integer in base-64
  notation -- i.e. as an array of 6 six-bit chunks   */
  
  long i, j, k, sum;
  longer mult, newseed;
  double x;

  mult[0] = 13;   /* these four statements set the multiplier */
  mult[1] = 24;   /* -- they are its "digits" in a base-64    */
  mult[2] = 22;   /*    notation: 1664525 = 6*64^3+22*64^2    */
  mult[3] = 6;    /*                         +24*64+13        */
  for (i = 0; i <= 5; i++)
    newseed[i] = 0;
  for (i = 0; i <= 5; i++) {  /* do the multiplication piecewise */
    sum = newseed[i];
    k = i;
    if (i > 3)
      k = 3;
    for (j = 0; j <= k; j++)
      sum += mult[j] * seed[i - j];
    newseed[i] = sum;
    for (j = i; j <= 4; j++) {
      newseed[j + 1] += newseed[j] / 64;
      newseed[j] &= 63;
    }
  }
  memcpy(seed, newseed, sizeof(longer));        /* new seed replaces old one */
  seed[5] &= 3;          /* from the new seed, get a floating point fraction */
  x = 0.0;
  for (i = 0; i <= 5; i++)
    x = x / 64.0 + seed[i];
  x /= 4.0;
  return x;
}  /* randum */


void randumize(longer seed, long *enterorder)
{ /* randomize input order of species -- randomly permute array enterorder */
  long i, j, k;

  for (i = 0; i < spp; i++) {
    j = (long)(randum(seed) * (i+1));
    k = enterorder[j];
    enterorder[j] = enterorder[i];
    enterorder[i] = k;
  }
} /* randumize */


double normrand(longer seed)
{/* standardized Normal random variate */
  double x;

  x = randum(seed)+randum(seed)+randum(seed)+randum(seed)
       + randum(seed)+randum(seed)+randum(seed)+randum(seed)
       + randum(seed)+randum(seed)+randum(seed)+randum(seed)-6.0;
  return(x);
} /* normrand */ 

/************* User configuration **************/

void initseed(long *inseed, long *inseed0, longer seed)
{ /* input random number seed */
  long i, loopcount;

  loopcount = 0;
  do {
    printf("Random number seed (must be odd)?\n");
    scanf("%ld%*[^\n]", inseed);
    getchar();
    countup(&loopcount, 10);
  } while (((*inseed) < 0) || ((*inseed) & 1) == 0);
  *inseed0 = *inseed;
  for (i = 0; i <= 5; i++)
    seed[i] = 0;
  i = 0;
  do {
    seed[i] = *inseed & 63;
    *inseed /= 64;
    i++;
  } while (*inseed != 0);
}  /*initseed*/


void initjumble(long *inseed, long *inseed0, longer seed, long *njumble)
{ /* input number of jumblings for jumble option */
  long loopcount;

  initseed(inseed, inseed0, seed);
  loopcount = 0;
  do {
    printf("Number of times to jumble?\n");
    fflush(stdout);
    scanf("%ld%*[^\n]", njumble);
    getchar();
    countup(&loopcount, 10);
  } while ((*njumble) < 1);
}  /*initjumble*/


void initoutgroup(long *outgrno, long spp)
{ /* input outgroup number */
  long loopcount;
  boolean done;

  loopcount = 0;
  do {
    printf("Type number of the outgroup:\n");
    fflush(stdout);
    scanf("%ld%*[^\n]", outgrno);
    getchar();
    done = (*outgrno >= 1 && *outgrno <= spp);
    if (!done) {
      printf("BAD OUTGROUP NUMBER: %ld\n", *outgrno);
      printf("  Must be in range 1 - %ld\n", spp);
    }
    countup(&loopcount, 10);
  } while (done != true);
}  /*initoutgroup*/


void initthreshold(double *threshold)
{ /* input threshold for threshold parsimony option */
  long loopcount;
  boolean done;

  loopcount = 0;
  do {
    printf("What will be the threshold value?\n");
    fflush(stdout);
    scanf("%lf%*[^\n]", threshold);
    getchar();
    done = (*threshold >= 1.0);
    if (!done)
      printf("BAD THRESHOLD VALUE:  it must be greater than 1\n");
    else
      *threshold = (long)(*threshold * 10.0 + 0.5) / 10.0;
    countup(&loopcount, 10);
  } while (done != true);
}  /*initthreshold*/


void initcatn(long *categs)
{ /* initialize category number for rate categories */
  long loopcount;

  loopcount = 0;
  *categs = 0;
  do {
    printf("Number of categories (1-%d)?\n", maxcategs);
    fflush(stdout);
    scanf("%ld%*[^\n]", categs);
    getchar();
    countup(&loopcount, 10);
  } while (*categs > maxcategs || *categs < 1);
}  /*initcatn*/


void initcategs(long categs, double *rate)
{ /* initialize category rates for HMM rates */
  long i, loopcount, scanned;
  char line[100], rest[100];
  boolean done;

  loopcount = 0;
  for (;;){
    printf("Rate for each category? (use a space to separate)\n");
    fflush(stdout);
    getstryng(line);
    done = true;
    for (i = 0; i < categs; i++){
      scanned = sscanf(line,"%lf %[^\n]", &rate[i],rest);
      if ((scanned < 2 && i < (categs - 1)) ||
       (scanned < 1 && i == (categs - 1))){
     printf("Please enter exactly %ld values.\n", categs);
     done = false;
     break;
      }
      strcpy(line, rest);
    }
    if (done)
      break;
    countup(&loopcount, 100);
  }
}  /*initcategs*/


void initprobcat(long categs, double *probsum, double *probcat)
{ /* input probabilities of rate categores for HMM rates */
  long i, loopcount, scanned;
  boolean done;
  char line[100], rest[100];

  loopcount = 0;
  do {
    printf("Probability for each category?");
    printf(" (use a space to separate)\n");
    fflush(stdout);
    getstryng(line);
    done = true;
    for (i = 0; i < categs; i++){
      scanned = sscanf(line, "%lf %[^\n]", &probcat[i], rest);
      if ((scanned < 2 && i < (categs - 1)) ||
       (scanned < 1 && i == (categs - 1))){
     done = false;
     printf("Please enter exactly %ld values.\n", categs);
     break;}
      strcpy(line, rest);
    }
    if (!done)
      continue;
    *probsum = 0.0;
    for (i = 0; i < categs; i++)
      *probsum += probcat[i];
    if (fabs(1.0 - (*probsum)) > 0.001) {
      done = false;
      printf("Probabilities must add up to");
      printf(" 1.0, plus or minus 0.001.\n");
    }
    countup(&loopcount, 100);
  } while (!done);
}  /*initprobcat*/

/************ Math utility functions ********/
void lgr(long m, double b, raterootarray lgroot)
{ /* For use by initgammacat.  Get roots of m-th Generalized Laguerre
     polynomial, given roots of (m-1)-th, these are to be
     stored in lgroot[m][] */
  long i;
  double upper, lower, x, y;
  boolean dwn;   /* is function declining in this interval? */

  if (m == 1) {
    lgroot[1][1] = 1.0+b;
  } else {
    dwn = true;
    for (i=1; i<=m; i++) {
      if (i < m) {
        if (i == 1)
          lower = 0.0;
        else
          lower = lgroot[m-1][i-1];
        upper = lgroot[m-1][i];
      } else {                 /* i == m, must search above */
        lower = lgroot[m-1][i-1];
        x = lgroot[m-1][m-1];
        do {
          x = 2.0*x;
          y = glaguerre(m, b, x);
        } while ((dwn && (y > 0.0)) || ((!dwn) && (y < 0.0)));
        upper = x;
      }
      while (upper-lower > 0.000000001) {
        x = (upper+lower)/2.0;
        if (glaguerre(m, b, x) > 0.0) {
          if (dwn)
            lower = x;
          else
            upper = x;
        } else {
          if (dwn)
            upper = x;
          else
            lower = x;
        }        
      }
      lgroot[m][i] = (lower+upper)/2.0;
      dwn = !dwn;                /* switch for next one */
    }
  }
} /* lgr */


double logfac (long n)
{ /* log(n!) values were calculated with Mathematica
     with a precision of 30 digits */
  long i;
  double x;

  switch (n)
    {
    case 0:
      return 0.;
    case 1:
      return 0.;
    case 2:
      return 0.693147180559945309417232121458;
    case 3:
      return 1.791759469228055000812477358381;
    case 4:
      return 3.1780538303479456196469416013;
    case 5:
      return 4.78749174278204599424770093452;
    case 6:
      return 6.5792512120101009950601782929;
    case 7:
      return 8.52516136106541430016553103635;
    case 8:
      return 10.60460290274525022841722740072;
    case 9:
      return 12.80182748008146961120771787457;
    case 10:
      return 15.10441257307551529522570932925;
    case 11:
      return 17.50230784587388583928765290722;
    case 12:
      return 19.98721449566188614951736238706;
    default:
      x = 19.98721449566188614951736238706;
      for (i = 13; i <= n; i++)
        x += log(i);
      return x;
    }
} /* logfac */

                        
double glaguerre(long m, double b, double x)
{ /* Generalized Laguerre polynomial computed recursively.
     For use by initgammacat */
  long i;
  double gln, glnm1, glnp1; /* L_n, L_(n-1), L_(n+1) */

  if (m == 0)
    return 1.0;
  else {
    if (m == 1)
      return 1.0 + b - x;
    else {
      gln = 1.0+b-x;
      glnm1 = 1.0;
      for (i=2; i <= m; i++) {
        glnp1 = ((2*(i-1)+b+1.0-x)*gln - (i-1+b)*glnm1)/i;
        glnm1 = gln;
        gln = glnp1;
      }
      return gln;
    }
  }
} /* glaguerre */


void initlaguerrecat(long categs, double alpha, double *rate, double *probcat)
{ /* calculate rates and probabilities to approximate Gamma distribution
     of rates with "categs" categories and shape parameter "alpha" using
     rates and weights from Generalized Laguerre quadrature */
  long i;
  raterootarray lgroot; /* roots of GLaguerre polynomials */
  double f, x, xi, y;

  alpha = alpha - 1.0;
  lgroot[1][1] = 1.0+alpha;
  for (i = 2; i <= categs; i++)
    lgr(i, alpha, lgroot);                   /* get roots for L^(a)_n */
  /* here get weights */
  /* Gamma weights are (1+a)(1+a/2) ... (1+a/n)*x_i/((n+1)^2 [L_{n+1}^a(x_i)]^2)  */
  f = 1;
  for (i = 1; i <= categs; i++)
    f *= (1.0+alpha/i);
  for (i = 1; i <= categs; i++) {
    xi = lgroot[categs][i];
    y = glaguerre(categs+1, alpha, xi);
    x = f*xi/((categs+1)*(categs+1)*y*y);
    rate[i-1] = xi/(1.0+alpha);
    probcat[i-1] = x;
  }
} /* initlaguerrecat */


double hermite(long n, double x)
{ /* calculates hermite polynomial with degree n and parameter x */
  /* seems to be unprecise for n>13 -> root finder does not converge*/
  double h1 = 1.;
  double h2 = 2. * x;
  double xx = 2. * x;
  long i;

  for (i = 1; i < n; i++) {
    xx = 2. * x * h2 - 2. * (i) * h1;
    h1 = h2;
    h2 = xx;
  }
  return xx;
} /* hermite */


void root_hermite(long n, double *hroot)
{ /* find roots of Hermite polynmials */
  long z;
  long ii;
  long start;

  if (n % 2 == 0) {
    start = n/2;
    z = 1;
  } else {
    start = n/2 + 1;
    z=2;
    hroot[start-1] = 0.0;
  }
  for (ii = start; ii < n; ii++) {         /* search only upwards*/
    hroot[ii] = halfroot(hermite, n, hroot[ii-1]+EPSILON, 1./n);
    hroot[start - z] = -hroot[ii];
    z++;
  }
} /* root_hermite */


double halfroot(double (*func)(long m, double x), long n, double startx,
                double delta)
{ /* searches from the bound (startx) only in one direction
     (by positive or negative delta, which results in
     other-bound=startx+delta)
     delta should be small.
     (*func) is a function with two arguments  */
  double xl;            /* lower x value */
  double xu;            /* upper x value */
  double xm = 0.0;
  double fu;
  double fl;
  double fm = 100000.;
  double gradient;
  boolean dwn = false;

  /* decide if we search above or below startx and escapes to trace back
     to the starting point that most often will be
     the root from the previous calculation */
  if (delta < 0) {
    xu = startx;
    xl = xu + delta;
  } else {
    xl = startx;
    xu = xl + delta;
  }
  delta = fabs(delta);
  fu = (*func)(n, xu);
  fl = (*func)(n, xl);
  gradient = (fl-fu)/(xl-xu);
  while(fabs(fm) > EPSILON) {        /* is root outside of our bracket?*/
    if ((fu<0.0 && fl<0.0) || (fu>0.0 && fl > 0.0)) {
      xu += delta;
      fu = (*func)(n, xu);
      fl = (*func)(n, xl);
      gradient = (fl-fu)/(xl-xu);
      dwn = (gradient < 0.0) ? true : false;
    } else {
      xm = xl - fl / gradient;
      fm = (*func)(n, xm);
      if (dwn) {
        if (fm > 0.) {
          xl = xm;
          fl = fm;
        } else {
          xu = xm;
          fu = fm;
        }
      } else {
        if (fm > 0.) {
          xu = xm;
          fu = fm;
        } else {
          xl = xm;
          fl = fm;
        }
      }
      gradient = (fl-fu)/(xl-xu);
    }
  }
  return xm;
} /* halfroot */


void hermite_weight(long n, double * hroot, double * weights)
{
  /* calculate the weights for the hermite polynomial at the roots
     using formula from Abramowitz and Stegun chapter 25.4.46 p.890 */
  long i;
  double hr2;
  double numerator;

  numerator = exp(0.6931471805599 * ( n-1.) + logfac(n)) / (n*n);
  for (i = 0; i < n; i++) {
    hr2 = hermite(n-1, hroot[i]);
    weights[i] = numerator / (hr2*hr2);
  }
} /* hermiteweight */


void inithermitcat(long categs, double alpha, double *rate, double *probcat)
{ /* calculates rates and probabilities */
  long i;
  double *hroot;
  double std;

  std = SQRT2 /sqrt(alpha);
  hroot = (double *) Malloc((categs+1) * sizeof(double));
  root_hermite(categs, hroot);         /* calculate roots */
  hermite_weight(categs, hroot, probcat);  /* set weights */
  for (i=0; i<categs; i++) {           /* set rates */
    rate[i] = 1.0 + std * hroot[i];
    probcat[i] = probcat[i];
  }
  free(hroot);
} /* inithermitcat */


void initgammacat (long categs, double alpha, double *rate, double *probcat)
{ /* calculate rates and probabilities to approximate Gamma distribution
   of rates with "categs" categories and shape parameter "alpha" using
   rates and weights from Generalized Laguerre quadrature or from
   Hermite quadrature */

  if (alpha >= 100.0)
    inithermitcat(categs, alpha, rate, probcat);
  else
    initlaguerrecat(categs, alpha, rate, probcat);
} /* initgammacat */


void inithowmany(long *howmanny, long howoften)
{/* input how many cycles */
  long loopcount;

  loopcount = 0;
  do { 
    printf("How many cycles of %4ld trees?\n", howoften);
    fflush(stdout);
    scanf("%ld%*[^\n]", howmanny);
    getchar();
    countup(&loopcount, 10);
  } while (*howmanny <= 0);
}  /*inithowmany*/


void inithowoften(long *howoften)
{ /* input how many trees per cycle */
  long loopcount;

  loopcount = 0;
  do {
    printf("How many trees per cycle?\n");
    fflush(stdout);
    scanf("%ld%*[^\n]", howoften);
    getchar();
    countup(&loopcount, 10);
  } while (*howoften <= 0);
}  /*inithowoften*/


void initlambda(double *lambda)
{ /* input patch length parameter for autocorrelated HMM rates */
  long loopcount;

  loopcount = 0;
  do {
    printf("Mean block length of sites having the same rate (greater than 1)?\n");
    fflush(stdout);
    scanf("%lf%*[^\n]", lambda);
    getchar();
    countup(&loopcount, 10);
  } while (*lambda <= 1.0);
  *lambda = 1.0 / *lambda;
}  /*initlambda*/


void initfreqs(double *freqap, double *freqcp, double *freqgp, double *freqtp)
{ /* input frequencies of the four bases */
  char *str;
  double freqa, freqc, freqg, freqt, sum;
  long loopcount;

  assert(freqap && freqcp && freqgp && freqtp);

  loopcount = 0;
  for (;;) {
    printf("Base frequencies for A, C, G, T/U ?\n");
    fflush(stdout);
    str = fgetline(stdin);
    if ( sscanf(str, "%lf%lf%lf%lf", &freqa, &freqc, &freqg, &freqt) == 4 ) {
      if (    freqa >= 0.0 && freqc >= 0.0 
           && freqg >= 0.0 && freqt >= 0.0 ) {
        sum = freqa + freqc + freqg + freqt;
        freqa /= sum;
        freqc /= sum;
        freqg /= sum;
        freqt /= sum;

        if (fabs(sum - 1.0) >= 1.0e-3) {
          printf("Normalized base frequencies are:\n%.3f %.3f %.3f %.3f\n"
              "(press enter)", freqa, freqc, freqg, freqt);
          fflush(stdout);
          fgetline(stdin);
        }
        break;        /* input OK */
      }
      else {
        printf("ERROR: Frequencies cannot be negative.\n");
      }
    }
    else {
      printf("ERROR: Please enter four numbers separated by spaces.\n");
    }
    countup(&loopcount, 10);
  }
  *freqap = freqa;
  *freqcp = freqc;
  *freqgp = freqg;
  *freqtp = freqt;
}  /* initfreqs */


void initratio(double *ttratio)
{ /* input transition/transversion ratio */
  long loopcount;

  loopcount = 0;
  do {
    printf("Transition/transversion ratio?\n");
    fflush(stdout);
    scanf("%lf%*[^\n]", ttratio);
    getchar();
    countup(&loopcount, 10);
  } while (*ttratio < 0.0);
}  /* initratio */


void initpower(double *power)
{
  printf("New power?\n");
  fflush(stdout);
  scanf("%lf%*[^\n]", power);
  getchar();
}  /*initpower*/


void initdatasets(long *datasets)
{
  /* handle multi-data set option */
  long loopcount;
  boolean done;

  loopcount = 0;
  do {
    printf("How many data sets?\n");
    fflush(stdout);
    scanf("%ld%*[^\n]", datasets);
    getchar();
    done = (*datasets > 1);
      if (!done)
     printf("Bad data sets number:  it must be greater than 1\n");
    countup(&loopcount, 10);
  } while (!done);
} /* initdatasets */


void justweights(long *datasets)
{
  /* handle multi-data set option by weights */
  long loopcount;
  boolean done;

  loopcount = 0;
  do {
    printf("How many sets of weights?\n");
    fflush(stdout);
    scanf("%ld%*[^\n]", datasets);
    getchar();
    done = (*datasets >= 1);
      if (!done)
     printf("BAD NUMBER:  it must be greater than 1\n");
    countup(&loopcount, 10);
  } while (!done);
} /* justweights */


void initterminal(boolean *ibmpc, boolean *ansi)
{
  /* handle terminal option */
  if (*ibmpc) {
    *ibmpc = false;
    *ansi = true;
  } else if (*ansi)
      *ansi = false;
    else
      *ibmpc = true;
}  /*initterminal*/


void initnumlines(long *screenlines)
{
  long loopcount;

  loopcount = 0;
  do {
    *screenlines = readlong("Number of lines on screen?\n");
    countup(&loopcount, 10);
  } while (*screenlines <= 12);
}  /*initnumlines*/


void newline(FILE *filename, long i, long j, long k)
{
  /* go to new line if i is a multiple of j, indent k spaces */
  long m;

  if ((i - 1) % j != 0 || i <= 1)
    return;
  putc('\n', filename);
  for (m = 1; m <= k; m++)
    putc(' ', filename);
}  /* newline */


/************* Sequence file routines **************/

void inputnumbers(long *spp, long *chars, long *nonodes, long n)
{
  /* Read numbers of species and characters from first line of a data set.
   * Return the results in *spp and *chars, respectively. Also returns
   * (*spp * 2 - n)  in *nonodes */

  if (fscanf(infile, "%ld%ld", spp, chars) != 2 || *spp <= 0 || *chars <= 0) {
    printf(
    "ERROR: Unable to read the number of species or characters in data set\n");
    printf(
      "The input file is incorrect (perhaps it was not saved text only).\n");
  }
  *nonodes = *spp * 2 - n;
}  /* inputnumbers */


void inputnumbers2(long *spp, long *nonodes, long n)
{
  /* read species number */

  if (fscanf(infile, "%ld", spp) != 1 || *spp <= 0) {
    printf("ERROR: Unable to read the number of species in data set\n");
    printf(
      "The input file is incorrect (perhaps it was not saved text only).\n");
  }
  fprintf(outfile, "\n%4ld Populations\n", *spp);
  *nonodes = *spp * 2 - n;
}  /* inputnumbers2 */


void samenumsp(long *chars, long ith)
{
  /* check if spp is same as the first set in other data sets */
  long cursp, curchs;

  if (eoln(infile)) 
    scan_eoln(infile);
  fscanf(infile, "%ld%ld", &cursp, &curchs);
  if (cursp != spp) {
    printf(
         "\n\nERROR: Inconsistent number of species in data set %ld\n\n", ith);
    exxit(-1);
  }
  *chars = curchs;
} /* samenumsp */


void samenumsp2(long ith)
{
  /* check if spp is same as the first set in other data sets */
  long cursp;

  if (eoln(infile)) 
    scan_eoln(infile);
  if (fscanf(infile, "%ld", &cursp) != 1) {
    printf("\n\nERROR: Unable to read number of species in data set %ld\n",
      ith);
    printf(
      "The input file is incorrect (perhaps it was not saved text only).\n");
    exxit(-1);
  }
  if (cursp != spp) {
    printf(
      "\n\nERROR: Inconsistent number of species in data set %ld\n\n", ith);
    exxit(-1);
  }
} /* samenumsp2 */


void readoptions(long *extranum, const char *options)
{ /* read option characters from input file */
  Char ch;

  while (!(eoln(infile))) {
    ch = gettc(infile);
    uppercase(&ch);
    if (strchr(options, ch) != NULL)
     (* extranum)++;
    else if (!(ch == ' ' || ch == '\t')) {
      printf("BAD OPTION CHARACTER: %c\n", ch);
      exxit(-1);
    }
  }
  scan_eoln(infile);
}  /* readoptions */


void matchoptions(Char *ch, const char *options)
{  /* match option characters to those in auxiliary options line in restriction
    * site data file */

  *ch = gettc(infile);
  uppercase(ch);
  if (strchr(options, *ch) == NULL) {
    printf("ERROR: Incorrect auxiliary options line");
    printf(" which starts with %c\n", *ch);
    exxit(-1);
  }
}  /* matchoptions */


void headings(long chars, const char *letters1, const char *letters2)
{
  long i, j;

  putc('\n', outfile);
  j = nmlngth + (chars + (chars - 1) / 10) / 2 - 5;
  if (j < nmlngth - 1)
    j = nmlngth - 1;
  if (j > 37)
    j = 37;
  fprintf(outfile, "Name");
  for (i = 1; i <= j; i++)
    putc(' ', outfile);
  fprintf(outfile, "%s\n", letters1);
  fprintf(outfile, "----");
  for (i = 1; i <= j; i++)
    putc(' ', outfile);
  fprintf(outfile, "%s\n\n", letters2);
}  /* headings */


void initname(long i)
{
  /* read in species name */
  boolean gotatab;
  long j;

  gotatab = false;
  for (j = 0; j < nmlngth; j++) {
    if (eoff(infile) || eoln(infile)){
      printf("\n\nERROR: end-of-line or end-of-file");
      printf(" in the middle of species name for species %ld\n\n", i+1);
      exxit(-1);
    }
    if (!gotatab) {  /* if no tab character has been read yet */
      nayme[i][j] = gettc(infile);
      if ((nayme[i][j] == '(') || (nayme[i][j] == ')') || (nayme[i][j] == ':')
        || (nayme[i][j] == ',') || (nayme[i][j] == ';') || (nayme[i][j] == '[')
        || (nayme[i][j] == ']')) {
    printf("\nERROR: Species name may not contain characters ( ) : ; , [ ] \n");
    printf("       In the name of species number %ld at position number %ld\n",
         i+1, j+1);
       printf("       there is character %c\n\n", nayme[i][j]);
        exxit(-1);
      }
      if (nayme[i][j] == '\t') {  /* check to see if was a tab character */
        nayme[i][j] = ' ';
        gotatab = true;
      }
    }
    else {  /* once a tab character has been seen, blank-fill */
      nayme[i][j] = ' ';
    }
  }
} /* initname */


/*********** Weight file routines **********/

void inputweights(long chars, steptr weight, boolean *weights)
{
  /* input the character weights, 0-9 and A-Z for weights 0 - 35 */
  Char ch;
  long i;

  for (i = 0; i < chars; i++) {
    do {
      if (eoln(weightfile)) 
        scan_eoln(weightfile);
      ch = gettc(weightfile);
      if (ch == '\n')
        ch = ' ';
    } while (ch == ' ');
    weight[i] = 1;
    if (isdigit(ch))
      weight[i] = ch - '0';
    else if (isalpha(ch)) {
      uppercase(&ch);
      weight[i] = ch - 'A' + 10;
    } else {
      printf("\n\nERROR: Bad weight character: %c\n\n", ch);
      exxit(-1);
    }
  }
  scan_eoln(weightfile);
  *weights = true;
}  /* inputweights */


void inputweights2(long a, long b, long *weightsum,
        steptr weight, boolean *weights, const char *prog)
{
  /* input the character weights, 0 or 1 */
  Char ch;
  long i;

  *weightsum = 0;
  for (i = a; i < b; i++) {
    do {
      if (eoln(weightfile))
        scan_eoln(weightfile);
      ch = gettc(weightfile);
    } while (ch == ' ');
    weight[i] = 1;
    if (ch == '0' || ch == '1')
      weight[i] = ch - '0';
    else {
      printf("\n\nERROR: Bad weight character: %c -- ", ch);
      printf("weights in %s must be 0 or 1\n", prog);
      exxit(-1);
    }
    *weightsum += weight[i];
  }
  *weights = true;
  scan_eoln(weightfile);
}  /* inputweights2 */


void printweights(FILE *filename, long inc, long chars,
        steptr weight, const char *letters)
{
  /* print out the weights of sites */
  long i, j;
  boolean letterweights;

  letterweights = false;
  for (i = 0; i < chars; i++)
    if (weight[i] > 9)
      letterweights = true;
  fprintf(filename, "\n    %s are weighted as follows:", letters);
  if (letterweights)
    fprintf(filename, " (A = 10, B = 11, etc.)\n");
  else
    putc('\n', filename);
  for (i = 0; i < chars; i++) {
    if (i % 60 == 0) {
      putc('\n', filename);
      for (j = 1; j <= nmlngth + 3; j++)
        putc(' ', filename);
    }
    if (weight[i+inc] < 10)
      fprintf(filename, "%ld", weight[i + inc]);
    else
      fprintf(filename, "%c", 'A'-10+(int)weight[i + inc]);
    if ((i+1) % 5 == 0 && (i+1) % 60 != 0)
      putc(' ', filename);
  }
  fprintf(filename, "\n\n");
}  /* printweights */

/************* Category file routines ***************/

void inputcategs(long a, long b, steptr category, long categs, const char *prog)
{
  /* input the categories, 1-9 */
  Char ch;
  long i;

  for (i = a; i < b; i++) {
    do {
      if (eoln(catfile)) 
        scan_eoln(catfile);
      ch = gettc(catfile);
    } while (ch == ' ');
    if ((ch >= '1') && (ch <= ('0'+categs)))
      category[i] = ch - '0';
    else {
      printf("\n\nERROR: Bad category character: %c", ch);
      printf(" -- categories in %s are currently 1-%ld\n", prog, categs);
      exxit(-1);
    }
  }
  scan_eoln(catfile);
}  /* inputcategs */


void printcategs(FILE *filename, long chars, steptr category,
                 const char *letters)
{
  /* print out the sitewise categories */
  long i, j;

  fprintf(filename, "\n    %s are:\n", letters);
  for (i = 0; i < chars; i++) {
    if (i % 60 == 0) {
      putc('\n', filename);
      for (j = 1; j <= nmlngth + 3; j++)
        putc(' ', filename);
    }
    fprintf(filename, "%ld", category[i]);
    if ((i+1) % 10 == 0 && (i+1) % 60 != 0)
      putc(' ', filename);
  }
  fprintf(filename, "\n\n");
}  /* printcategs */

/*********** Factors file routines **********/

void inputfactors(long chars, Char *factor, boolean *factors)
{
  /* reads the factor symbols */
  long i;

  for (i = 0; i < (chars); i++) {
    if (eoln(factfile)) 
      scan_eoln(factfile);
    factor[i] = gettc(factfile);
    if (factor[i] == '\n')
      factor[i] = ' ';
  }
  scan_eoln(factfile);
  *factors = true;
}  /* inputfactors */


void printfactors(FILE *filename, long chars, Char *factor, const char *letters)
{
  /* print out list of factor symbols */
  long i;

  fprintf(filename, "Factors%s:\n\n", letters);
  for (i = 1; i <= nmlngth - 5; i++)
    putc(' ', filename);
  for (i = 1; i <= (chars); i++) {
    newline(filename, i, 55, nmlngth + 3);
    putc(factor[i - 1], filename);
    if (i % 5 == 0)
      putc(' ', filename);
  }
  putc('\n', filename);
}  /* printfactors */


void findtree(boolean *found, long *pos, long nextree, long *place, bestelm *bestrees)
{
  /* finds tree given by array place in array bestrees by binary search */
  /* used by dnacomp, dnapars, dollop, mix, & protpars */
  long i, lower, upper;
  boolean below, done;

  below = false;
  lower = 1;
  upper = nextree - 1;
  (*found) = false;
  while (!(*found) && lower <= upper) {
    (*pos) = (lower + upper) / 2;
    i = 3;
    done = false;
    while (!done) {
      done = (i > spp);
      if (!done)
        done = (place[i - 1] != bestrees[(*pos) - 1].btree[i - 1]);
      if (!done)
        i++;
    }
    (*found) = (i > spp);
    if (*found)
      break;
    below = (place[i - 1] <  bestrees[(*pos )- 1].btree[i - 1]);
    if (below)
      upper = (*pos) - 1;
    else
      lower = (*pos) + 1;
  }
  if (!(*found) && !below)
    (*pos)++;
}  /* findtree */


void addtree(long pos, long *nextree, boolean collapse, long *place, bestelm *bestrees)
{
  /* puts tree from array place in its proper position in array bestrees */
  /* used by dnacomp, dnapars, dollop, mix, & protpars */
  long i;

  for (i = *nextree - 1; i >= pos; i--){
    memcpy(bestrees[i].btree, bestrees[i - 1].btree, spp * sizeof(long));
    bestrees[i].gloreange = bestrees[i - 1].gloreange;
    bestrees[i - 1].gloreange = false;
    bestrees[i].locreange = bestrees[i - 1].locreange;
    bestrees[i - 1].locreange = false;
    bestrees[i].collapse = bestrees[i - 1].collapse;
  }
  for (i = 0; i < spp; i++)
    bestrees[pos - 1].btree[i] = place[i];
/*  bestrees[pos -1].gloreange = false;	*/
  bestrees[pos -1].collapse = false;
    
  (*nextree)++;
}  /* addtree */


long findunrearranged(bestelm *bestrees, long nextree, boolean glob)
{
  /* finds bestree with either global or local field false */
  long i;

  if (glob) {
    for (i = 0; i < nextree - 1; i++)
      if (!bestrees[i].gloreange)
        return i;
  } else {
    for (i = 0; i < nextree - 1; i++)
      if (!bestrees[i].locreange)
        return i;
  }
  return -1;
} /* findunrearranged */


void shellsort(double *a, long *b, long n)
{ /* Shell sort keeping a, b in same order */
  /* used by dnapenny, dolpenny, penny, and threshml */
  long gap, i, j, itemp;
  double rtemp;

  gap = n / 2;
  while (gap > 0) {
    for (i = gap + 1; i <= n; i++) {
      j = i - gap;
      while (j > 0) {
     if (a[j - 1] > a[j + gap - 1]) {
       rtemp = a[j - 1];
       a[j - 1] = a[j + gap - 1];
       a[j + gap - 1] = rtemp;
       itemp = b[j - 1];
       b[j - 1] = b[j + gap - 1];
       b[j + gap - 1] = itemp;
     }
     j -= gap;
      }
    }
    gap /= 2;
  }
}  /* shellsort */

/******** User trees ****************/

void getch(Char *c, long *parens, FILE *treefile)
{ /* get next nonblank character */

  do {
    if (eoln(treefile)) 
      scan_eoln(treefile);
    (*c) = gettc(treefile);

    if ((*c) == '\n' || (*c) == '\t')
      (*c) = ' ';
  } while ( *c == ' ' && !eoff(treefile) );
  if ((*c) == '(')
    (*parens)++;
  if ((*c) == ')')
    (*parens)--;
}  /* getch */


void findch(Char c, Char *ch, long which)
{ /* scan forward until find character c */
  boolean done;
  long dummy_parens;
  done = false;
  while (!done) {
    if (c == ',') {
      if (*ch == '(' || *ch == ')' || *ch == ';') {
        printf(
      "\n\nERROR in user tree %ld: unmatched parenthesis or missing comma\n\n",
          which);
        exxit(-1);
      } else if (*ch == ',')
        done = true;
    } else if (c == ')') {
      if (*ch == '(' || *ch == ',' || *ch == ';') {
        printf("\n\nERROR in user tree %ld: ", which);
        printf("unmatched parenthesis or non-bifurcated node\n\n");
        exxit(-1);
      } else {
        if (*ch == ')')
          done = true;
      }
    } else if (c == ';') {
      if (*ch != ';') {
        printf("\n\nERROR in user tree %ld: ", which);
        printf("unmatched parenthesis or missing semicolon\n\n");
        exxit(-1);
      } else
        done = true;
    }
    if (*ch != ')' && done)
      continue;
   getch(ch, &dummy_parens, intree);
  }
}  /* findch */


void processlength(double *valyew, double *divisor, Char *ch, 
        boolean *lengthIsNegative, FILE *treefile, long *parens)
{ /* read a branch length from a treefile */
  long digit, ordzero, exponent, exponentIsNegative;
  boolean pointread, hasExponent;

  ordzero = '0';
  *lengthIsNegative = false;
  pointread = false;
  hasExponent = false;
  exponentIsNegative = -1; // 3 states:  -1 = unassigned, 1 = true, 0 = false
  exponent = 0;
  *valyew = 0.0;
  *divisor = 1.0;
  getch(ch, parens, treefile);
  if ('+' == *ch)
    getch(ch, parens, treefile); // ignore leading '+', because "+1.2345" == "1.2345"
  else if ('-' == *ch)
    {
      *lengthIsNegative = true;
      getch(ch, parens, treefile);
    }
  digit = (long)(*ch - ordzero);
  while ( ((digit <= 9) && (digit >= 0)) || '.' == *ch || '-' == *ch
	  || '+' == *ch || 'E' == *ch || 'e' == *ch) {
    if ('.' == *ch)
      {
	if (!pointread)
	  pointread = true;
	else
	  {
	    printf("\n\nERROR: Branch length found with more than one \'.\' in it.\n\n");
	    exxit(-1);
	  }
      }
    else if ('+' == *ch)
      {
	if (hasExponent && -1 == exponentIsNegative)
	  exponentIsNegative = 0; // 3 states:  -1 = unassigned, 1 = true, 0 = false
	else
	  {
	    printf("\n\nERROR: Branch length found with \'+\' in an unexpected place.\n\n");
	    exxit(-1);
	  }
      }
    else if ('-' == *ch)
      {
	if (hasExponent && -1 == exponentIsNegative)
	  exponentIsNegative = 1; // 3 states:  -1 = unassigned, 1 = true, 0 = false
	else
	  {
	    printf("\n\nERROR: Branch length found with \'-\' in an unexpected place.\n\n");
	    exxit(-1);
	  }
      }
    else if ('E' == *ch || 'e' == *ch)
      {
	if (!hasExponent)
	  hasExponent = true;
	else
	  {
	    printf("\n\nERROR: Branch length found with more than one \'E\' in it.\n\n");
	    exxit(-1);
	  }
      }
    else {
      if (!hasExponent)
	{
	  *valyew = *valyew * 10.0 + digit;
	  if (pointread)
	    *divisor *= 10.0;
	}
      else
	exponent = 10*exponent + digit;
    }
    getch(ch, parens, treefile);
    digit = (long)(*ch - ordzero);
  }
  if (hasExponent)
    {
      if (exponentIsNegative)
	*divisor *= pow(10.,(double)exponent);
      else
	*divisor /= pow(10.,(double)exponent);
    }
  if (*lengthIsNegative)
    *valyew = -(*valyew);
}  /* processlength */


void commentskipper(FILE *intree, long *bracket)
{ /* skip over comment bracket contents in reading tree */
  char c;
  
  c = gettc(intree);
  
  while (c != ']') {
    
    if(feof(intree)) {
      printf("\n\nERROR: Unmatched comment brackets\n\n");
      exxit(-1);
    }

    if(c == '[') {
      (*bracket)++;
      commentskipper(intree, bracket);
    }
    c = gettc(intree);
  }
  (*bracket)--;
}  /* commentskipper */


long countcomma(FILE *treefile, long *comma)
{
  /* Modified by Dan F. 11/10/96 */ 

  /* countcomma rewritten so it passes back both lparen+comma to allocate nodep
    and a pointer to the comma variable.  This allows the tree to know how many
    species exist, and the tips to be placed in the front of the nodep array */
  /* The next line inserted so this function leaves the file pointing
     to where it found it, not just re-winding it. */
  long orig_position = ftell(treefile);

  Char c;
  long  lparen = 0;
  long bracket = 0;
  (*comma) = 0;


  for (;;){
    c = getc(treefile);
    if (feof(treefile))
      break;
    if (c == ';')
      break;
    if (c == ',')
      (*comma)++;
    if (c == '(')
         lparen++;
    if (c == '[') {
      bracket++;
      commentskipper(treefile, &bracket);
    }
  }

  /* Don't just rewind, */
  /* rewind (treefile); */
  /* Re-set to where it pointed when the function was called */

  fseek (treefile, orig_position, SEEK_SET);

  return lparen + (*comma);
}  /*countcomma*/


long countsemic(FILE *treefile)
{ /* Used to determine the number of user trees.  Return
     either a: the number of semicolons in the file outside comments
     or b: the first integer in the file */
  Char c;
  long return_val, semic = 0;
  long bracket = 0;
  
  /* Eat all whitespace */
  c = gettc(treefile);
  while ((c == ' ')  ||
      (c == '\t') ||
      (c == '\n')) {
    c = gettc(treefile);
  }

  /* Then figure out if the first non-white character is a digit; if
     so, return it */
  if (isdigit (c)) {
    ungetc(c, treefile);
    fscanf(treefile, "%ld", &return_val);
  } else {

    /* Loop past all characters, count the number of semicolons
       outside of comments */
    for (;;){
      c = fgetc(treefile);
      if (feof(treefile))
     break;
      if (c == ';')
     semic++;
      if (c == '[') {
     bracket++;
     commentskipper(treefile, &bracket);
      }
    }
    return_val = semic;
  }

  rewind (treefile);
  return return_val;
}  /* countsemic */


/************* Memory management *************/

void memerror()
{
  printf("Error allocating memory\n");
  exxit(-1);
}  /* memerror */


void odd_malloc(long x)
{ /* error message if attempt to malloc too little or too much memory */
  printf ("ERROR: a function asked for an inappropriate amount of memory:");
  printf ("  %ld bytes\n", x);
  printf ("       This can mean one of two things:\n");
  printf ("       1.  The input file is incorrect");
  printf (" (perhaps it was not saved as Text Only),\n");
  printf ("       2.  There is a bug in the program.\n");
  printf ("       Please check your input file carefully.\n");
  printf ("       If it seems to be a bug, please mail joe@gs.washington.edu\n");
  printf ("       with the name of the program, your computer system type,\n");
  printf ("       a full description of the problem, and with the input data file.\n");

  /* abort() can be used to crash */
  
  exxit(-1);
}


MALLOCRETURN *mymalloc(long x)
{ /* wrapper for malloc, allowing error message if too little, too much */
  MALLOCRETURN *new_block;

  if ((x <= 0) ||
      (x > TOO_MUCH_MEMORY))
    odd_malloc(x);

  new_block = (MALLOCRETURN *)calloc(1, x);

  if (!new_block) {
    memerror();
    return (MALLOCRETURN *) new_block;
  } else
    return (MALLOCRETURN *) new_block;
} /* mymalloc */


void hookup(node *p, node *q)
{ /* hook together two nodes */

  /* IMPORTANT -- does not change branch lengths. Other routines
   * expect them to be as they were, and update them later */
  p->back = q;
  q->back = p;
}  /* hookup */


void link_trees(long local_nextnum, long nodenum, long local_nodenum,
        pointarray nodep)
{
  if(local_nextnum == 0)
    hookup(nodep[nodenum], nodep[local_nodenum]);
  else if(local_nextnum == 1)
      hookup(nodep[nodenum], nodep[local_nodenum]->next);
    else if(local_nextnum == 2)
        hookup(nodep[nodenum], nodep[local_nodenum]->next->next);
      else
        printf("Error in Link_Trees()");
} /* link_trees() */


void allocate_nodep(pointarray *nodep, FILE *treefile, long  *precalc_tips)  
{ /* pre-compute space and allocate memory for nodep */

  long numnodes;      /* returns number commas & (    */
  long numcom = 0;        /* returns number commas */
  
  numnodes = countcomma(treefile, &numcom) + 1;
  *nodep      = (pointarray)Malloc(2*numnodes*sizeof(node *));

  (*precalc_tips) = numcom + 1;        /* this will be used in placing the
                                          tip nodes in the front region of
                                          nodep.  Used for species check?  */
} /* allocate_nodep -plc */



long take_name_from_tree (Char *ch, Char *str, FILE *treefile)
{
  /* This loop reads a name from treefile and stores it in *str.
     Returns the length of the name string. str must be at
     least MAXNCH bytes, but no effort is made to null-terminate
     the string. Underscores and newlines are converted to spaces.
     Characters beyond MAXNCH are discarded. */

  long name_length = 0;

  do {
    if ((*ch) == '_')
      (*ch) = ' ';
    if ( name_length < MAXNCH )
      str[name_length++] = (*ch);
    if (eoln(treefile)) 
      scan_eoln(treefile);
    (*ch) = gettc(treefile);
    if (*ch == '\n')
      *ch = ' ';
  } while ( strchr(":,)[;", *ch) == NULL );

  return name_length;
}  /* take_name_from_tree */


void match_names_to_data (Char *str, pointarray treenode, node **p, long spp)
{
  /* This loop matches names taken from treefile to indexed names in
     the data file */

  boolean found;
  long i, n;

  n = 1;  
  do {
    found = true;
    for (i = 0; i < nmlngth; i++) {
      found = (found && ((str[i] == nayme[n - 1][i]) ||
        (((nayme[n - 1][i] == '_') && (str[i] == ' ')) ||
        ((nayme[n - 1][i] == ' ') && (str[i] == '\0')))));
    }
    
    if (found)
      *p = treenode[n - 1];
    else
      n++;

  } while (!(n > spp || found));
  
  if (n > spp) {
    printf("\n\nERROR: Cannot find species: ");
    for (i = 0; (str[i] != '\0') && (i < MAXNCH); i++)
      putchar(str[i]);
    printf(" in data file\n\n");
    exxit(-1);
  }
}  /* match_names_to_data */


void addelement(node **p, node *q, Char *ch, long *parens, FILE *treefile,
        pointarray treenode, boolean *goteof, boolean *first, pointarray nodep,
        long *nextnode, long *ntips, boolean *haslengths, node **grbg,
        initptr initnode, boolean unifok, long maxnodes)
{
  /* Recursive procedure adds nodes to user-defined tree
     This is the main (new) tree-reading procedure */

  node *pfirst;
  long i, len = 0, nodei = 0;
  boolean notlast;
  Char str[MAXNCH+1];
  node *r;
  long furs = 0;

  if ((*ch) == '(') {
    (*nextnode)++;          /* get ready to use new interior node */
    nodei = *nextnode;      /* do what needs to be done at bottom */
    if ( maxnodes != -1 && nodei > maxnodes) {
      printf("ERROR in input tree file: Attempting to allocate too\n"); 
      printf("many nodes. This is usually caused by a unifurcation.\n");  
      printf("To use this tree with this program  use Retree to read\n");
      printf("and write this tree.\n");
      exxit(-1);
    }
    
    /* do what needs to be done at bottom */
    (*initnode)(p, grbg, q, len, nodei, ntips,
                  parens, bottom, treenode, nodep, str, ch, treefile);
    pfirst      = (*p);
    notlast = true;
    while (notlast) {          /* loop through immediate descendants */
      furs++;
      (*initnode)(&(*p)->next, grbg, q,
                   len, nodei, ntips, parens, nonbottom, treenode,
                   nodep, str, ch, treefile);
                       /* ... doing what is done before each */
      r = (*p)->next;
      getch(ch, parens, treefile);      /* look for next character */
      
       /* handle blank names */
      if((*ch) == ',' || (*ch) == ':'){
        ungetc((*ch), treefile);
        *ch = 0;
      } else if((*ch)==')'){
        ungetc((*ch), treefile);
        (*parens)++;
        *ch = 0;
      }
 
      addelement(&(*p)->next->back, (*p)->next, ch, parens, treefile,
        treenode, goteof, first, nodep, nextnode, ntips,
        haslengths, grbg, initnode, unifok, maxnodes);  

      (*initnode)(&r, grbg, q, len, nodei, ntips,
                    parens, hslength, treenode, nodep, str, ch, treefile);
                           /* do what is done after each about length */
      *p = r;                         /* make r point back to p */

      if ((*ch) == ')') {
        notlast = false;
        do {
          getch(ch, parens, treefile);
        } while ((*ch) != ',' && (*ch) != ')' &&
           (*ch) != '[' && (*ch) != ';' && (*ch) != ':');
      }
    }
    if ( furs <= 1 && !unifok ) {
      printf("ERROR in input tree file: A Unifurcation was detetected.\n");
      printf("To use this tree with this program use retree to read and");
      printf(" write this tree\n");
      exxit(-1);
    }
    
    (*p)->next = pfirst;
    (*p)       = pfirst;

  } else if ((*ch) != ')') {       /* if it's a species name */
    for (i = 0; i < MAXNCH+1; i++)   /* fill string with nulls */
      str[i] = '\0';

    len = take_name_from_tree (ch, str, treefile);  /* get the name */

    if ((*ch) == ')')
      (*parens)--;         /* decrement count of open parentheses */
    (*initnode)(p, grbg, q, len, nodei, ntips,
                  parens, tip, treenode, nodep, str, ch, treefile);
                              /* do what needs to be done at a tip */
  } else
    getch(ch, parens, treefile);
  if (q != NULL)
    hookup(q, (*p));                    /* now hook up */
  (*initnode)(p, grbg, q, len, nodei, ntips, 
                parens, iter, treenode, nodep, str, ch, treefile);
                              /* do what needs to be done to variable iter */
  if ((*ch) == ':')
    (*initnode)(p, grbg, q, len, nodei, ntips, 
                  parens, length, treenode, nodep, str, ch, treefile);
                                   /* do what needs to be done with length */
  else if ((*ch) != ';' && (*ch) != '[')
    (*initnode)(p, grbg, q, len, nodei, ntips, 
                  parens, hsnolength, treenode, nodep, str, ch, treefile);
                             /* ... or what needs to be done when no length */
  if ((*ch) == '[')
    (*initnode)(p, grbg, q, len, nodei, ntips,
                  parens, treewt, treenode, nodep, str, ch, treefile);
                             /* ... for processing a tree weight */
  else if ((*ch) == ';')     /* ... and at end of tree */
    (*initnode)(p, grbg, q, len, nodei, ntips,
                  parens, unittrwt, treenode, nodep, str, ch, treefile);
}  /* addelement */


void treeread (FILE *treefile, node **root, pointarray treenode,
        boolean *goteof, boolean *first, pointarray nodep, 
        long *nextnode, boolean *haslengths, node **grbg, initptr initnode,
        boolean unifok, long maxnodes)
{
  /* read in user-defined tree and set it up */
  /* Eats blank lines and everything up to the first open paren, then
   * calls the recursive function addelement, which builds the
   * tree and calls back to initnode. */
  char  ch;
  long parens = 0;
  long ntips = 0;
  
  (*goteof) = false;
  (*nextnode) = spp;

  /* eat blank lines */
  while (eoln(treefile) && !eoff(treefile)) 
    scan_eoln(treefile);

  if (eoff(treefile)) {
    (*goteof) = true;
    return;
  } 

  getch(&ch, &parens, treefile);

  while (ch != '(') {
    /* Eat everything in the file (i.e. digits, tabs) until you
       encounter an open-paren */
    getch(&ch, &parens, treefile);
  }
  (*haslengths) = true; 
  addelement(root, NULL, &ch, &parens, treefile,
         treenode, goteof, first, nodep, nextnode, &ntips,
         haslengths, grbg, initnode, unifok, maxnodes);

  /* Eat blank lines and end of current line*/
  do {
    scan_eoln(treefile);
  }
  while (eoln(treefile) && !eoff(treefile));

  (*first) = false;
  if (parens != 0) {
    printf("\n\nERROR in tree file: unmatched parentheses\n\n");
    exxit(-1);
  }
}  /* treeread */


void addelement2(node *q, Char *ch, long *parens, FILE *treefile,
        pointarray treenode, boolean lngths, double *trweight, boolean *goteof,
        long *nextnode, long *ntips, long no_species, boolean *haslengths,
        boolean unifok, long maxnodes)
{
  /* recursive procedure adds nodes to user-defined tree
     -- old-style bifurcating-only version */
  node *pfirst = NULL, *p;
  long i, len, current_loop_index;
  boolean notlast, minusread;
  Char str[MAXNCH];
  double valyew, divisor;
  long furs = 0; 

  if ((*ch) == '(') {

    current_loop_index = (*nextnode) + spp;
    (*nextnode)++;

    if ( maxnodes != -1 && current_loop_index > maxnodes) {
      printf("ERROR in intree file: Attempting to allocate too many nodes\n");
      printf("This is usually caused by a unifurcation.  To use this\n");
      printf("intree with this program  use retree to read and write\n"); 
      printf("this tree.\n");
      exxit(-1);
    }
    /* This is an assignment of an interior node */
    p = treenode[current_loop_index];
    pfirst = p;
    notlast = true;
    while (notlast) {
      furs++;
      /* This while loop goes through a circle (triad for
      bifurcations) of nodes */
      p = p->next;
      /* added to ensure that non base nodes in loops have indices */
      p->index = current_loop_index + 1;
      
      getch(ch, parens, treefile);
      
      addelement2(p, ch, parens, treefile, treenode, lngths, trweight,
        goteof, nextnode, ntips, no_species, haslengths, unifok, maxnodes);

      if ((*ch) == ')') {
        notlast = false;
        do {
          getch(ch, parens, treefile);
        } while ((*ch) != ',' && (*ch) != ')' &&
           (*ch) != '[' && (*ch) != ';' && (*ch) != ':');
      }
    }
    if ( furs <= 1 && !unifok ) {
      printf("ERROR in intree file: A Unifurcation was detected.\n");
      printf("To use this intree with this program use retree to read and");
      printf(" write this tree\n");
      exxit(-1);
    }

  } else if ((*ch) != ')') {
    for (i = 0; i < MAXNCH; i++) 
      str[i] = '\0';
    len = take_name_from_tree (ch, str, treefile);
    match_names_to_data (str, treenode, &p, spp);
    pfirst = p;
    if ((*ch) == ')')
      (*parens)--;
    (*ntips)++;
    strncpy (p->nayme, str, len);
  } else
    getch(ch, parens, treefile);
  
  if ((*ch) == '[') {    /* getting tree weight from last comment field */
    if (!eoln(treefile)) {
      fscanf(treefile, "%lf", trweight);
      getch(ch, parens, treefile);
      if (*ch != ']') {
        printf("\n\nERROR: Missing right square bracket\n\n");
        exxit(-1);
      }
      else {
        getch(ch, parens, treefile);
        if (*ch != ';') {
          printf("\n\nERROR: Missing semicolon after square brackets\n\n");
          exxit(-1);
        }
      }
    }
  }
  else if ((*ch) == ';') {
    (*trweight) = 1.0 ;
    if (!eoln(treefile))
      printf("WARNING: tree weight set to 1.0\n");
  }
  else
    (*haslengths) = ((*haslengths) && q == NULL);
  
  if (q != NULL)
    hookup(q, pfirst);
  /*if (q != NULL) {
    if (q->branchnum < pfirst->branchnum)
      pfirst->branchnum = q->branchnum;
    else
      q->branchnum = pfirst->branchnum;
  } FIXME check if we need this for restml */

  if ((*ch) == ':') {
    processlength(&valyew, &divisor, ch,
       &minusread, treefile, parens);
    if (q != NULL) {
      if (!minusread)
        q->oldlen = valyew / divisor;
      else
        q->oldlen = 0.0;
      if (lngths) {
        q->v = valyew / divisor;
        q->back->v = q->v;
        q->iter = false;
        q->back->iter = false;
        q->back->iter = false;
    }
    }
  }
  
}  /* addelement2 */


void treeread2 (FILE *treefile, node **root, pointarray treenode,
        boolean lngths, double *trweight, boolean *goteof,
        boolean *haslengths, long *no_species, boolean unifok, long maxnodes)
{
  /* read in user-defined tree and set it up
     -- old-style bifurcating-only version */
  char  ch;
  long parens = 0;
  long ntips = 0;
  long nextnode;
  
  (*goteof) = false;
  nextnode = 0;

  /* Eats all blank lines at start of file */
  while (eoln(treefile) && !eoff(treefile)) 
    scan_eoln(treefile);

  if (eoff(treefile)) {
    (*goteof) = true;
    return;
  } 

  getch(&ch, &parens, treefile);

  while (ch != '(') {
    /* Eat everything in the file (i.e. digits, tabs) until you
       encounter an open-paren */
    getch(&ch, &parens, treefile);
  }

  addelement2(NULL, &ch, &parens, treefile, treenode, lngths, trweight,
      goteof, &nextnode, &ntips, (*no_species), haslengths, unifok, maxnodes);
  (*root) = treenode[*no_species];

  /*eat blank lines */
  while (eoln(treefile) && !eoff(treefile)) 
    scan_eoln(treefile);

  (*root)->oldlen = 0.0;

  if (parens != 0) {
    printf("\n\nERROR in tree file:  unmatched parentheses\n\n");
    exxit(-1);
  }
}  /* treeread2 */


void exxit(int exitcode)
{ /* Terminate the program with exit code exitcode.
   * On Windows, supplying a nonzero exitcode will print a message and wait
   * for the user to hit enter. */

#if defined(WIN32) || defined(OSX_CARBON)
  if (exitcode == 0)
    exit (exitcode);
  else {
    puts ("Hit Enter or Return to close program.");
    fgetline(stdin);
#endif
#ifdef WIN32
    phyRestoreConsoleAttributes();
#endif
#if defined(WIN32) || defined(OSX_CARBON)
    exit (exitcode);
  }
#else
  exit(exitcode);
#endif

} /* exxit */


char gettc(FILE* file) 
{ /* Return the next character in file.
   * If EOF is reached, print an error and die.
   * DOS ('\r\n') and Mac ('\r') newlines are returned as a single '\n'. */
  int ch;

  ch=getc(file);

  if ( ch == EOF )
    EOF_error();

  if ( ch == '\r' ) {
    ch = getc(file);
    if ( ch != '\n' )
      ungetc(ch, file);
    ch = '\n';
  }
  return ch;
} /* gettc */


/************* More tree functions **********/

void unroot(tree *t, long nonodes) 
{
  /* used by fitch, restml and contml */
  if (t->root->back == NULL) { 
    if (t->root->next->back->tip)
      t->root = t->root->next->next->back;
    else  t->root = t->root->next->back;
  }
  if (t->root->next->back == NULL) {
    if (t->root->back->tip)
      t->root = t->root->next->next->back;
    else t->root = t->root->back;
  }
  if (t->root->next->next->back == NULL)  {
    if (t->root->back->tip)
      t->root = t->root->next->back;
    else t->root = t->root->back;
  }
    

  unroot_r(t->root, t->nodep, nonodes);
  unroot_r(t->root->back, t->nodep, nonodes);
}


void unroot_here(node* root, node** nodep, long nonodes)
{
  node* tmpnode;
  double newl;
  /* used by unroot */
  /* assumes bifurcation this is ok in the programs that use it */
 
  
  newl = root->next->oldlen + root->next->next->oldlen;
  root->next->back->oldlen = newl;
  root->next->next->back->oldlen = newl;

  newl = root->next->v + root->next->next->v;
  root->next->back->v = newl;
  root->next->next->back->v = newl;

  root->next->back->back=root->next->next->back;
  root->next->next->back->back = root->next->back;
  
  while ( root->index != nonodes ) {
    tmpnode = nodep[ root->index ];
    nodep[root->index] = root;
    root->index++;
    root->next->index++;
    root->next->next->index++;
    nodep[root->index - 2] = tmpnode;
    tmpnode->index--;
    tmpnode->next->index--;
    tmpnode->next->next->index--;
  }
}


void unroot_r(node* p, node** nodep, long nonodes) 
{
  /* used by unroot */
  node *q;
  
  if ( p->tip) return;

  q = p->next;
  while ( q != p ) {
    if (q->back == NULL)
      unroot_here(q, nodep, nonodes);
    else unroot_r(q->back, nodep, nonodes);
    q = q->next;
  }
}


void destruct_tree(tree* t) 
{ /* returns a tree such that there are no branches, and the free fork nodes
     go on the stacks */
  long j, nsibs;
  node *q, *p;

  while ( !Slist_isempty(t->free_forks) ) t->get_fork(t);
  for ( j = 0; j < t->nonodes ; j++ )  {
    p = t->nodep[j];
    p->back = NULL;
    /* EWFIX.REMOVE
    p->initialized = false;
    */
    if ( j < spp ) continue;
    
    /* Get rid of extra furcations if necessary */
    for ( nsibs = count_sibs(p); nsibs > 2; nsibs-- ) {
      q = p->next->next;
      t->release_forknode(t, p->next);
      p->next = q;
    }
    
    p->initialized = false;
    p->next->initialized = false;
    p->next->next->initialized = false;
    p->back = NULL;
    p->next->back = NULL;
    p->next->next->back = NULL;
    
    t->release_fork(t, p); 
  }

}


void generic_tree_free(tree *t)
{
  long i;
  node *p,*q,*r;


  free(t->free_forks);  /* EWFIX -- doesn't get everything on list */
  free(t->free_fork_nodes);

  for ( i = 0 ; i < NLRSAVES ; i++ )
    t->lrsaves[i]->free(&(t->lrsaves[i]));
  free(t->lrsaves);
  t->temp_p->free(&(t->temp_p));
  t->temp_q->free(&(t->temp_q));

  for ( i = 0 ; i < t->nonodes ; i++ ) {
    p = t->nodep[i];
    if ( i >= spp ) {
      q = p->next;
      while ( q != p ) {
        r = q->next;
        q->free(&q);
        q = r;
      }
    }
    p->free(&t->nodep[i]);
  }
  free(t->nodep);
  free(t);
}


void rooted_tree_init(tree* t, long nonodes, long spp)
{/* a few extra things for a rooted tree*/
  generic_tree_init(t, nonodes, spp);
  t->globrearrange = rooted_globrearrange;
  t->insert_ = rooted_tree_insert_;
  t->re_move = rooted_tree_re_move;
  t->locrearrange = rooted_locrearrange;
  t->save_lr_nodes = rooted_tree_save_lr_nodes;
  t->restore_lr_nodes = rooted_tree_restore_lr_nodes;
}


void generic_tree_init(tree* t, long nonodes, long spp)
{
  long i, j;
  node *q,*p;

  /* these functions may be customized for each program */
  if ( t->release_fork == NULL )
    t->release_fork = generic_tree_release_fork;
  if ( t->get_fork == NULL )
    t->get_fork = generic_tree_get_fork;
  if ( t->release_forknode == NULL )
    t->release_forknode = generic_tree_release_forknode;

  t->spp = spp;
  t->nonodes = nonodes;
  t->nodep = Malloc(nonodes* sizeof(node *));
  for ( i = 0 ; i < spp ; i++ ) {
    t->nodep[i] = functions.new_node(true, i+1);
  }
  for ( i = spp ; i <  nonodes ; i++ ) {
    q = NULL;
    for ( j = 1 ; j <= 3 ; j++ ) {
      p = functions.new_node(false, i + 1 );
      p->next = q;
      q = p;
    }
    p->next->next->next = p;
    t->nodep[i] = p;
  }

  /* Create garbage lists */
  t->free_forks = Slist_new();
  t->free_fork_nodes = Slist_new();

  /* Put all nodes on garbage lists by "releasing" them */
  for ( i = nonodes - 1 ; i >= spp ; i-- ) {
    t->release_fork(t, t->nodep[i]);
  }

  t->lrsaves = Malloc(NLRSAVES * sizeof(node*));
  for ( i = 0 ; i < NLRSAVES ; i++ )
    t->lrsaves[i] = functions.new_node(false,0);
  t->temp_p = functions.new_node(false,0);
  t->temp_q = functions.new_node(false,0);

  t->addtraverse = generic_tree_addtraverse;
  t->globrearrange = generic_globrearrange;
  t->free = generic_tree_free;
  t->copy = generic_tree_copy;
  t->smoothall = (tree_smoothall_t)no_op;
  t->root = t->nodep[0];
  t->root = NULL;
  t->score = UNDEFINED;
  t->locrearrange = generic_unrooted_locrearrange;
  t->save_lr_nodes = unrooted_tree_save_lr_nodes;
  t->restore_lr_nodes = unrooted_tree_restore_lr_nodes; 
  t->save_traverses = generic_tree_save_traverses;
  t->restore_traverses = generic_tree_restore_traverses;
  t->nuview = generic_tree_nuview;
  t->evaluate = generic_tree_evaluate;
  t->insert_ = generic_tree_insert_;
  t->get_forknode = generic_tree_get_forknode;
  t->re_move = generic_tree_re_move;
  t->try_insert_ = generic_tree_try_insert_;
  t->tree_print_f = generic_tree_print;
  t->do_branchl_on_insert_f = generic_do_branchl_on_insert;
  t->do_branchl_on_re_move_f = generic_do_branchl_on_re_move;

  t->tree_good_f = generic_tree_good;   /* debug  ?? */
  t->node_good_f = generic_node_good;
  t->fork_good_f = generic_fork_good;
} /* generic_tree-init */


tree* generic_tree_new(long nonodes, long spp) 
{
  tree* t = Malloc(sizeof(tree));
  generic_tree_init(t, nonodes, spp); 
  return t;
}

void generic_tree_print(tree * t)
{
    long nodeIndex;
    printf("-----------------------------------------------\n");
    printf("tree %p spp = %ld ; nonodes = %ld root = %p\n",
        t,t->spp,t->nonodes,t->root);
    for(nodeIndex=0;nodeIndex<t->nonodes;nodeIndex++)
    {
        node * p = t->nodep[nodeIndex];
        printf("---\nnodep[%ld]: %p",nodeIndex,p);
        if(p->tip) printf(" TIP");
        if(p == t->root) printf(" ROOT");
        printf("\n");
        if(p)
        {
            p->fork_print_f(p);
        }
    }
}

boolean generic_tree_good(tree *t)
{
    long nodeIndex;
    for(nodeIndex = 0; nodeIndex < t->nonodes; nodeIndex++)
    {
        node * n = t->nodep[nodeIndex];
        if( n->tip )
        {
            boolean thisNodeGood = t->node_good_f(t,n);
            if (!thisNodeGood) return false;
        }
        else
        {
            boolean thisNodeGood = t->fork_good_f(t,n);
            if (!thisNodeGood) return false;
        }
    }
    return true;
}

boolean generic_fork_good(tree *t, node * n)
{
    boolean firstTime = true;
    boolean hasNullBack = false;
    boolean hasGoodBack = false;
    node * p = n;
    while ( firstTime || (p != n)) 
    {
        if(p == NULL)
        {
            assert(false);
            return false;
        }
        else
        {
            if (p->back == NULL) 
                hasNullBack = true;
            else
                hasGoodBack = true;
            boolean nodeGood = t->node_good_f(t,p);
            if ( !nodeGood )
            {
                assert(false);
                return false;
            }
            p = p->next;
        }
        firstTime = false;
    }
    return true;
}


boolean generic_node_good(tree *t, node * n)
{
    if ( n->back != NULL)
    {
        boolean edgesEqual = (n->back->v == n->v);
        assert(edgesEqual);
        if ( !edgesEqual) return false;
    }
    return true;
}


void rooted_globrearrange(tree* curtree, boolean progress, boolean thorough)
{
  /* does global rearrangements */
  tree *globtree, *oldtree, *priortree, *bestree;
  int i;
  node *where,*sib_ptr,*qwhere;
  double oldbestyet;
  int success = false;
  boolean succeeded = true;
  double bestyet;
  boolean multf;
 
  /* FIXME should do the "Doing global rearrangements" printf here instead of 
   * outside of this function in every program */

  bestree = functions.tree_new(curtree->nonodes, curtree->spp); 
  globtree = functions.tree_new(curtree->nonodes, curtree->spp);
  priortree = functions.tree_new(curtree->nonodes, curtree->spp);
  oldtree = functions.tree_new(curtree->nonodes, curtree->spp);

  while ( succeeded ) {
    if (progress) {
      printf("   ");
      fflush(stdout);
    }

    succeeded = false;
    bestyet = oldbestyet = curtree->score;
    curtree->copy(curtree, globtree);
    curtree->copy(curtree, oldtree);
   
    for ( i = 0 ; i < curtree->nonodes ; i++ ) {
      bestyet = curtree->score;
      sib_ptr  = curtree->nodep[i];

      if ( progress && ((i - spp) % (( curtree->nonodes / 72 ) + 1 ) == 0 ))
        putchar('.');
      fflush(stdout);

      if (sib_ptr->index == curtree->root->index)
        continue;
      if ( sib_ptr->back == NULL ) /* this implies unused node */
        continue; /* probably because of multifurcation */
    
      curtree->re_move(curtree, sib_ptr,&where, true); 
      curtree->copy(curtree, priortree); 
      qwhere = where;
        
      succeeded  = curtree->addtraverse(curtree, sib_ptr, curtree->root, true, 
          &qwhere, &bestyet, bestree, priortree, thorough, &multf); 
      if ( thorough ) {
        if ( where != qwhere && bestyet > globtree->score) {
            bestree->copy(bestree, globtree);
            success = true;
        }
      } else {
        if ( succeeded && where != qwhere) {
          curtree->insert_(curtree, sib_ptr, qwhere, true, multf);
          curtree->smoothall(curtree, where);
          success = true;
          curtree->copy(curtree, globtree);
        }
      }
      oldtree->copy(oldtree, curtree);
      oldtree->copy(oldtree, bestree);
    }
    globtree->copy(globtree, curtree);
    globtree->copy(globtree, bestree);
    succeeded = success && globtree->score > oldbestyet;
    
    if (progress)
      putchar('\n');
  }
  
  bestree->free(bestree);
  priortree->free(priortree);
  globtree->free(globtree);
  oldtree->free(oldtree);
}


void generic_globrearrange(tree* curtree, boolean progress, boolean thorough)
{ /* does global rearrangements */
  tree *globtree, *oldtree, *priortree, *bestree;
  int i, j, k, num_sibs, num_sibs2;
  node *where,*sib_ptr,*sib_ptr2, *qwhere;
  double oldbestyet, bestyet;
  int success = false;
  boolean succeeded = true;
  boolean multf;
  node* removed;
 
  if ( progress )  {
    printf("Doing global rearrangements\n");
    printf("  !");
    for ( i = 0 ; i < curtree->nonodes ; i++)
      printf("-");
    printf("!\n");
        phyFillScreenColor();
	fflush(stdout);
  }

  bestree = functions.tree_new(curtree->nonodes, curtree->spp); 
  globtree = functions.tree_new(curtree->nonodes, curtree->spp);
  priortree = functions.tree_new(curtree->nonodes, curtree->spp);
  oldtree = functions.tree_new(curtree->nonodes, curtree->spp);
  
  while ( succeeded ) {
    succeeded = false;
    curtree->smoothall(curtree, curtree->root);
    bestyet = oldbestyet = curtree->score;

    if (progress) {
      printf("   ");
      fflush(stdout);
    }

    curtree->copy(curtree, globtree);
    curtree->copy(curtree, oldtree);
    
    for ( i = 0 ; i < curtree->nonodes ; i++ ) {
      sib_ptr  = curtree->nodep[i];
      if ( sib_ptr->tip )
        num_sibs = 0;
      else
        num_sibs = count_sibs(sib_ptr);

      if ( progress && ((i - spp) % (( curtree->nonodes / 72 ) + 1 ) == 0 ))
        putchar('.');
      fflush(stdout);

      for ( j = 0 ; j <= num_sibs ; j++ ) {
        sib_ptr = curtree->nodep[i];
        for ( k = 0 ; k < j ; k++ )
          sib_ptr = sib_ptr->next;
        if ( sib_ptr->back == NULL || sib_ptr->back->tip ) 
          continue;

        removed = sib_ptr;
          

        curtree->re_move(curtree, removed,&where, true); 
        curtree->smoothall(curtree, where);
        curtree->copy(curtree, priortree); 
        qwhere = where;
       
        if ( where->tip) {
          num_sibs2 = 0;
          sib_ptr2 = where->back;
        }
        else {
          num_sibs2 = count_sibs(where);
          sib_ptr2 = where;
        }
        for ( k = 0 ; k <= num_sibs2 ; k++ ) {
          succeeded = curtree->addtraverse(curtree, removed, sib_ptr2->back,
              true, &qwhere,&bestyet, bestree, priortree, thorough,&multf) 
              || succeeded;
          sib_ptr2 = sib_ptr2->next;
        }
        if ( !thorough) {
            if (succeeded && qwhere != where && qwhere != where->back && 
                bestyet > oldbestyet) {
              curtree->insert_(curtree, removed, qwhere, true, multf);
              curtree->smoothall(curtree, where);
              success = true;
              curtree->copy(curtree, globtree);
            }
            else
              globtree->copy(globtree, curtree);
            
        }
          else {
            if ( qwhere && where != qwhere && bestyet > globtree->score) {
              bestree->copy(bestree, globtree);
              success = true;
            }
            oldtree->copy(oldtree, curtree);
            oldtree->copy(oldtree, bestree);
        }
      }
    }
    globtree->copy(globtree, curtree);
    globtree->copy(globtree, bestree);
    globtree->copy(globtree, oldtree);
    succeeded = success && globtree->score > oldbestyet;
    
    if (progress)
      putchar('\n');
  }
  
  bestree->free(bestree);
  priortree->free(priortree);
  globtree->free(globtree);
  oldtree->free(oldtree);
}


boolean generic_tree_addtraverse(tree* t, node *p, node*q, boolean contin, 
    node **qwherein, double* bestyet, tree* bestree, tree* priortree,
    boolean thorough, boolean* multf)
{ /* try adding p at q, proceed recursively through tree */
  node *sib_ptr;
  boolean succeeded= false;

  succeeded = t->try_insert_(t, p, q, qwherein, bestyet, bestree, priortree,
      thorough, multf);

  if (!q->tip && contin) {
    for ( sib_ptr = q->next ; q != sib_ptr ; sib_ptr = sib_ptr->next) {
      succeeded = generic_tree_addtraverse(t, p, sib_ptr->back, 
          contin, qwherein, bestyet, bestree, priortree, thorough, multf) 
        || succeeded;
    }
  }
  if (contin && q == t->root && t->root->back && 
      (t->root->back->tip == false)) { 
    /* we need to go both ways, in an unrooted tree */
    for ( sib_ptr = t->root->back->next ; t->root->back != sib_ptr ; 
        sib_ptr = sib_ptr->next) {
      succeeded = generic_tree_addtraverse(t, p, sib_ptr->back, 
          contin, qwherein, bestyet, bestree, priortree, thorough, multf) 
        || succeeded;
    }
  }
  return succeeded;
}


#ifdef WIN32
void phySaveConsoleAttributes()
{
  if ( GetConsoleScreenBufferInfo(hConsoleOutput, &savecsbi) )
    savecsbi_valid = true;
} /* PhySaveConsoleAttributes */


void phySetConsoleAttributes()
{
  hConsoleOutput = GetStdHandle(STD_OUTPUT_HANDLE);

  if ( hConsoleOutput == INVALID_HANDLE_VALUE )
    hConsoleOutput == NULL;

  if ( hConsoleOutput != NULL ) {
    phySaveConsoleAttributes();

    SetConsoleTextAttribute(hConsoleOutput, 
      BACKGROUND_GREEN | BACKGROUND_BLUE | BACKGROUND_INTENSITY);
  }
} /* phySetConsoleAttributes */


void phyRestoreConsoleAttributes()
{
  COORD coordScreen = { 0, 0 };
  DWORD cCharsWritten;
  DWORD dwConSize; 

  if ( savecsbi_valid ) {
    dwConSize = savecsbi.dwSize.X * savecsbi.dwSize.Y;

    SetConsoleTextAttribute(hConsoleOutput, savecsbi.wAttributes);

    FillConsoleOutputAttribute( hConsoleOutput, savecsbi.wAttributes,
           dwConSize, coordScreen, &cCharsWritten );
  }
} /* phyRestoreConsoleAttributes */


void phyFillScreenColor()
{
  COORD coordScreen = { 0, 0 };
  DWORD cCharsWritten;
  CONSOLE_SCREEN_BUFFER_INFO csbi; /* to get buffer info */ 
  DWORD dwConSize; 

  if ( GetConsoleScreenBufferInfo( hConsoleOutput, &csbi ) ) {
    dwConSize = csbi.dwSize.X * csbi.dwSize.Y;

    FillConsoleOutputAttribute( hConsoleOutput, csbi.wAttributes,
           dwConSize, coordScreen, &cCharsWritten );
  }
} /* phyFillScreenColor */


void phyClearScreen()
{
   COORD coordScreen = { 0, 0 };    /* here's where we'll home the
                                       cursor */ 
   DWORD cCharsWritten;
   CONSOLE_SCREEN_BUFFER_INFO csbi; /* to get buffer info */ 
   DWORD dwConSize;                 /* number of character cells in
                                       the current buffer */ 

   /* get the number of character cells in the current buffer */ 

   if ( GetConsoleScreenBufferInfo(hConsoleOutput, &csbi) ) {
       dwConSize = csbi.dwSize.X * csbi.dwSize.Y;

       /* fill the entire screen with blanks */ 

       FillConsoleOutputCharacter( hConsoleOutput, (TCHAR) ' ',
          dwConSize, coordScreen, &cCharsWritten );

       /* get the current text attribute */ 
       GetConsoleScreenBufferInfo( hConsoleOutput, &csbi );

       /* now set the buffer's attributes accordingly */ 
       FillConsoleOutputAttribute( hConsoleOutput, csbi.wAttributes,
          dwConSize, coordScreen, &cCharsWritten );

       /* put the cursor at (0, 0) */ 
       SetConsoleCursorPosition( hConsoleOutput, coordScreen );
  } 
} /* phyClearScreen */

#endif /* WIN32 */


void unrooted_tree_save_lr_nodes(tree* t, node* p, node* r)
{
  r->back->copy(r->back, t->lrsaves[0]);
  r->back->next->copy(r->back->next, t->lrsaves[1]);
  r->back->next->next->copy(r->back->next->next, t->lrsaves[2]);
  p->next->copy(p->next, t->lrsaves[3]);
  p->next->next->copy(p->next->next, t->lrsaves[4]);
  t->rb = r->back;
  t->rnb = r->back->next->back;
  t->rnnb = r->back->next->next->back;
}


void unrooted_tree_restore_lr_nodes(tree* t, node* p, node* r) 
{

  t->lrsaves[0]->copy(t->lrsaves[0], t->rb);
  t->lrsaves[1]->copy(t->lrsaves[1], t->rnb->back);
  t->lrsaves[2]->copy(t->lrsaves[2], t->rnnb->back);
  t->lrsaves[3]->copy(t->lrsaves[3], p->next);
  t->lrsaves[4]->copy(t->lrsaves[4], p->next->next);


  t->rb->back->v = t->rb->v;
  t->rnb->back->v = t->rnb->v;
  t->rnnb->back->v = t->rnnb->v;
  p->next->back->v = p->next->v;
  p->next->next->back->v = p->next->next->v;

  /* EWFIX.INIT */
  inittrav(t->rb);
  inittrav(t->rnb);
  inittrav(t->rnnb);
  inittrav(p->next);
  inittrav(p->next->next);

}


void generic_unrooted_locrearrange(tree* t, node* start, boolean thorough,
    tree* priortree, tree* bestree)
{ /* generic form of local rearrangement */
  double bestyet = t->evaluate(t, start, 0);
  boolean succeeded = true;
  while(succeeded) {
    succeeded = unrooted_tree_locrearrange_recurs(t, start->back, start, 
        &bestyet, thorough, priortree, bestree);
  }
} /* generic_unrooted_locrearrange */


boolean unrooted_tree_locrearrange_recurs(tree* t, node *p, node*pp, 
    double* bestyet, boolean thorough, tree* priortree, tree* bestree) 
{
  /* rearranges the tree locally moving pp around near p  */
  /* this function doesn't handle multifurcations */
  node *q, *r, *qwhere;
  boolean succeeded = false;
  boolean multf = false;
  double oldbestyet;
  qwhere = NULL;

  if (!p->tip && !p->back->tip)
  {
    /* Why are we setting t->score here? */
    /*    oldbestyet = t->score = *bestyet;  */

    /* Do we want to use the current tree instead? */
    /* oldbestyet = t->evaluate(t, t->root, 0);
     *
     * ...or use the former value of bestyet?  -IDR */
    oldbestyet = *bestyet;

    if (p->back->next != pp)
      r = p->back->next->back;
    else
      r = p->back->next->next->back;

            /* EWFIX.TREECHECK
            printf("EWFIX: before attempt\n");
            printf(" p : %p ; pp : %p ; r : %p\n",p,pp,r);
            t->tree_print_f(t);
            assert(t->tree_good_f(t));
            */

    if (!thorough) 
      t->save_lr_nodes(t, p, r);
    else
      t->copy(t, bestree);
    t->re_move(t, r, &q, false);

            /* EWFIX.TREECHECK
            printf("EWFIX: after re_move\n");
            t->tree_print_f(t);
            assert(t->tree_good_f(t));
            */

    if (thorough)
      t->copy(t, priortree);
    else
      qwhere = q;

    t->addtraverse(t, r, p->next, false,&qwhere, bestyet, bestree,
        priortree, thorough,&multf); 

            /* EWFIX.TREECHECK
            printf("EWFIX: after addtraverse\n");
            t->tree_print_f(t);
            assert(t->tree_good_f(t));
            */

    if(qwhere == q)
        // don't continue if we've already got a better tree
    {

        t->addtraverse(t, r, p->next->next, false,&qwhere, bestyet, bestree,
            priortree, thorough,&multf); 

            /* EWFIX.TREECHECK
            printf("EWFIX: after second addtraverse\n");
            t->tree_print_f(t);
            assert(t->tree_good_f(t));
            */
        }

    if (thorough)
      bestree->copy(bestree, t);
    else {
      if (qwhere == q ) {
        assert(*bestyet <= oldbestyet);
        t->insert_(t, r, qwhere, true, multf);

            /* EWFIX.TREECHECK
            printf("EWFIX: after re-insert\n");
            t->tree_print_f(t);
            assert(t->tree_good_f(t));
            */

        t->restore_lr_nodes(t, p, r);
        t->score = *bestyet;

            /* EWFIX.TREECHECK
            printf("EWFIX: after restore\n");
            t->tree_print_f(t);
            assert(t->tree_good_f(t));
            */
          }
      else {
        assert(*bestyet > oldbestyet);
        succeeded = true;
        t->insert_(t, r, qwhere, true, multf);

            /* EWFIX.TREECHECK
            printf("EWFIX: after insert of %p at %p\n",r,qwhere);
            t->tree_print_f(t);
            assert(t->tree_good_f(t));
            */

        t->smoothall(t, r->back);
        *bestyet = t->evaluate(t, p,0);

            /* EWFIX.TREECHECK
            printf("EWFIX: after smooth and evaluate\n");
            t->tree_print_f(t);
            assert(t->tree_good_f(t));
            */

/* debug        double otherbest = *bestyet;      JF:  is this needed? */
        assert(*bestyet = t->evaluate(t,t->root,0));
        /* EWFIX -- how different do we expect these two to be? */

      }
    }
    assert(oldbestyet <= *bestyet );
  }
  
  /* If rearrangements failed here, try subtrees, but stop when we find
   * one that improves the score. */
  if ( !p->tip && !succeeded ) {
    if ( unrooted_tree_locrearrange_recurs(t, p->next->back, p, bestyet,
        thorough, priortree, bestree) )
      return true;
    if ( unrooted_tree_locrearrange_recurs(t, p->next->next->back, p,
        bestyet, thorough, priortree, bestree) )
      return true;
  }
  return succeeded;
}
/* unrooted_tree_locrearrange_recurs */


/* generic_tree_save_traverses
 *
 * Saves the branch lengths for p and q (args to insert_) in t
 * This way, we can insert a fork above q and still recover
 * the original tree.
 */
void generic_tree_save_traverses(tree* t, node * p, node* q)
{
  p->copy(p,t->temp_p);
  q->copy(q,t->temp_q);
}


/* generic_tree_restore_traverses
 *
 * Restores branch legths to p and q (args to re_move) from
 * temp_p and temp_q nodes in t
 */
void generic_tree_restore_traverses(tree* t, node *p, node* q)
{
  t->temp_p->copy(t->temp_p,p);
  t->temp_q->copy(t->temp_q,q);
  inittrav(p);
  inittrav(q);
  if ( p->back )
  {
      p->back->v = p->v;
      inittrav(p->back);
  }
  if ( q->back )
  {
      q->back->v = q->v;
      inittrav(q->back);
  }
  /* EWFIX.INIT -- might be more correct to do all inittravs after ->v updates */
  /* EWFIX.TREECHECK
  printf("EWFIX restoring %p and %p\n\t",p,q);
  p->node_print_f(p);
  printf("\n\t");
  q->node_print_f(q);
  printf("\n");
  */


  // note: removed code to restore back links and release
  // fork node. this is now done as part of re_move
}


static void rooted_tryrearr(tree*t, node *p, boolean *success)
{
  /* evaluates one rearrangement of the tree.
    if the new tree has greater score than the old
    one sets success = TRUE and keeps the new tree.
    otherwise, restores the old tree */
  /* TODO relatively trivial to add a thorough mode */
  node *whereto, *forknode, *where;
  double oldlike, like; 

  p = t->nodep[p->index - 1];
  if (p == t->root)
    return;
  forknode = t->nodep[p->back->index - 1];
  if (forknode == t->root)
    return;
  oldlike = t->score;
   
  whereto = t->nodep[forknode->back->index - 1];
  t->save_lr_nodes(t, p, whereto);
  t->re_move(t, p, &where, false);
  t->insert_(t, p, whereto, false, false);
  like = t->evaluate(t, p, false);
  if (like - oldlike < LIKE_EPSILON)  {
    t->restore_lr_nodes(t, p, whereto);
    t->score = oldlike;
  } else {
    (*success) = true;
    t->smoothall(t, t->root);
  }
}  /* tryrearr */


static void rooted_repreorder(tree* t, node *p, boolean *success)
{
  /* traverses a binary tree, calling function mlk_tryrearr
    at a node before calling mlk_tryrearr at its descendants */
  node* q;
  if (p == NULL)
    return;
  rooted_tryrearr(t, p, success);
  if (p->tip)
    return;
  for ( q = p->next ; q != p && !(*success) ; q = q->next ) 
    rooted_repreorder(t, q->back, success);
}  /* repreorder */


void rooted_locrearrange(tree* t, node* start, boolean thorough,
    tree* priortree, tree* bestree)
{
  /*traverses the tree (preorder), finding any local
    rearrangement which increases the score.
    if traversal succeeds in increasing the tree's
    score, function rearrange runs traversal again  */

  boolean success;

  t->evaluate(t, start, 0); /* need to start of with a valid t->score */
  success = true;
  while (success) {
    success = false;
    rooted_repreorder(t, start, &success);
  }
}  /* rearrange */


void rooted_tree_save_lr_nodes(tree* t, node* p, node* whereto)
{
  node* forknode = t->nodep[p->back->index - 1];
  
  p->back->copy(p->back, t->lrsaves[0]);
  whereto->copy(whereto, t->lrsaves[1]);

  t->rnb = forknode->back;
  if ( p == forknode->next->back ) {
    t->onleft = false;
    t->rnnb = forknode->next->next->back;
  } else {
    t->onleft = true;
    t->rnnb = forknode->next->back;
  }
  whereto->initialized = false;
  p->back->initialized = false;
}


void rooted_tree_restore_lr_nodes(tree* t, node* p, node* whereto) 
{
  node* forknode = t->nodep[p->back->index - 1];

  if ( p == forknode->next->back )  {
    if (forknode->back != NULL)
      hookup( forknode->back, forknode->next->next->back);
    else {
      forknode->next->next->back->back = NULL;
      t->root = forknode->next->next->back;
    }
  } else {
    if ( forknode->back != NULL) 
      hookup( forknode->back, forknode->next->back);
    else {
      forknode->next->back->back = NULL;
      t->root = forknode->next->back;
    }
  }

  hookup(forknode, t->rnb);
  if ( t->onleft ) {
    hookup(forknode->next->next, p);
    hookup(forknode->next, t->rnnb);
  } else  {
    hookup(forknode->next, p);
    hookup(forknode->next->next, t->rnnb);
  }

  t->lrsaves[0]->copy(t->lrsaves[0], p->back);
  t->lrsaves[1]->copy(t->lrsaves[1], whereto);
}


void* pop(stack** oldstack)
{
  void* retval;
  stack* newstack;

  retval = (*oldstack)->data;
  newstack = (*oldstack)->next;
  free(*oldstack);
  *oldstack = newstack;
  return retval;
}


stack* push(stack* oldstack, void* newdata)
{
  stack* newstack;

  newstack = Malloc(sizeof(stack));
  newstack->data = newdata;
  newstack->next = oldstack;
  return newstack;
}


/* generic_tree_get_fork
 *
 * Pop a fork (ring of 3 nodes) off the free_forks stack, set initialized to
 * false on all, and return.
 */
node* generic_tree_get_fork(tree* t)
{
  node* retval;
  retval = Slist_pop(t->free_forks);
  retval->initialized = false;
  retval->next->initialized = false;
  retval->next->next->initialized = false;
  return retval;
}


void generic_tree_release_fork(tree* t, node* n)
{
  node *p;
  long sibs;

  /* we will only change the n pointer if we really need to */
  sibs = count_sibs(n);
  if ( sibs > 2 ) n = t->nodep[n->index  - 1];

  /* sibs initialized in the previous line */
  /* release forknodes until the fork is made of three forknodes */
  for ( ; sibs > 2 ; sibs-- ) {
    p = n->next;
    n->next = n->next->next;
    t->release_forknode(t, p);
  }


  /* remove any potentially contaminated data from remaining fork nodes */
  boolean firstTime = true;
  for ( p = n ; p != n || firstTime; p = p->next )
  {
     p->reinit(p);
     firstTime = false;
  }

  assert(count_sibs(n) == 2);
  Slist_push(t->free_forks, n);
}


void generic_tree_nuview(tree* t, node*p )
{
  node* sib_ptr;

  /* Recursive calls, should be called for all children */
  for ( sib_ptr = p->next ; sib_ptr != p ; sib_ptr = sib_ptr->next ) {
    if ( sib_ptr->back && !sib_ptr->back->tip && !sib_ptr->back->initialized)
    {
      t->nuview (t, sib_ptr->back);
    }
    else
    {
    }
  }
}


/* generic_tree_evaluate
 *
 * Updates nuviews for p and p->back in preparation for evaluation specific
 * to each program.
 */
double generic_tree_evaluate(tree *t, node* p, boolean dummy)
{
  node *q;

  q = p->back;
  if ( p->initialized == false && p->tip == false )
  {
      t->nuview((tree*)t, p);
  }
  if ( q  && q->initialized == false && q->tip == false ) 
  {
      t->nuview((tree*)t, q);
  }
  return 0;
}


void generic_tree_insert_(tree* t, node*p, node*q, boolean doinit, 
    boolean multf)
{ 
  node *newnode;

  if ( !multf ) {
    newnode = t->get_fork(t);

    assert(newnode->next->next->next==newnode);

    hookup(newnode, p);
    hookup(newnode->next->next, q->back);
    hookup(newnode->next, q);

    t->do_branchl_on_insert_f(t,newnode,q);

    assert( ! newnode->initialized );
    assert( ! newnode->next->initialized );
    assert( ! newnode->next->next->initialized );

    /* EWFIX.INIT
    if (doinit) {
    */
      inittrav(p);
      inittrav(p->back);
    /* EWFIX.INIT
    }
    */
  }
  else {
    newnode = t->get_forknode(t, q->index);
    newnode->next = q->next;
    q->next = newnode;
    hookup(newnode, p);      

    assert( ! newnode->initialized );

    if ( doinit ) {
      inittrav(p);
      inittrav(p->back);
    }
  }
}

/* split branch length when inserting */
void generic_do_branchl_on_insert(tree*t, node *fork, node* q) 
{
  /*
   * see ml.c for an example
   */
}


/* generic_tree_get_forknode
 *
 * Return an unused node with index i.
 *
 * If there are any nodes on the free_fork_nodes stack, one of these
 * is returned. Otherwise, create a new node and return it.
 */
node* generic_tree_get_forknode(tree* t, long i)
{
  node *p;
  if ( Slist_isempty(t->free_fork_nodes) )
    p = functions.new_node(0, i);
  else {
    p = Slist_pop(t->free_fork_nodes);
    p->init(p, 0, i);
  }
  return p;
}


void generic_tree_re_move(tree* t, node* item, node** where, boolean doinit)
{
  node *fork,*q,*p;
  long num_sibs;   
  
  fork = item->back;
  item->back = NULL;

  if ( item->tip && fork->tip ) {
    item->back = NULL;
    fork->back = NULL;
    return;
  }
  num_sibs = count_sibs(fork); 
  
  if ( num_sibs > 2 ) {
    for ( q = fork ; q->next != fork ; q = q->next)
      /* nothing */;
    
    q->next = fork->next;
    fork->next = NULL;
    fork->back = NULL;
    if ( t->nodep[fork->index - 1] == fork )
      t->nodep[fork->index - 1] = q;
    t->release_forknode(t, fork);
    if ( t->root == fork )
      t->root = q;
    if ( doinit ) {
      inittrav(q);
      for ( p = q->next ; p != q ; p = p->next )
        inittrav(p);
    }
    (*where) = q;
    
  } else {
    (*where) = fork->next->back;
    fork->next->back->back = fork->next->next->back;
    fork->next->next->back->back = fork->next->back;
    item->back = NULL;
    if ( fork == t->root )
      t->root = *where;
    if (t->root->tip ) t->root = t->root->back;

    t->do_branchl_on_re_move_f(t,item,*where);

    t->release_fork(t, fork);

    /* EWFIX.INIT -- might be in do_branchl_on_re_move ?? */
    if ( doinit ) {
      inittrav(*where);
      inittrav((*where)->back);
    }
  }
} /* generic_tree_re_move */

void generic_do_branchl_on_re_move(tree * t, node * p, node *q)
{
    /* see version in ml.c */
}


void generic_tree_release_forknode(tree* t, node* n)
{
  n->reinit(n);
  n->next = NULL;   // node_reinit(n) sets n->back to NULL
  Slist_push(t->free_fork_nodes, n);
} /* generic_tree_release_forknode */


boolean generic_tree_try_insert_(tree *t, node *p, node *q, node** qwherein, 
    double* bestyet, tree* bestree, tree* priortree, boolean thorough,
    boolean* multf)
{
  double like;
  boolean succeeded = false;
  node* dummy;

  t->insert_(t, p, q, true, false);
  like = t->evaluate(t, p, false);
  if (like > *bestyet + LIKE_EPSILON || *bestyet == UNDEFINED) {
    *bestyet = like;
    *qwherein = q;
    succeeded = true;
    *multf = false;
    if (thorough)
      t->copy(t, bestree);
  }
  if ( thorough )
    priortree->copy(priortree, t);
  else
    t->re_move(t, p, &dummy, false);
  
  return succeeded;
} /* generic_tree_try_insert_ */


void rooted_tree_insert_(tree* t, node* newtip, node* below, boolean doinit, 
    boolean multf)
  /* Insert node newtip into the tree above node below, adding a new fork
   * if necessary. If multf is TRUE, newtip is added as a new child of below,
   * without an additional fork.
   * 
   * TODO: implement the following:
   * If t->root is NULL, below is ignored, no fork is added, and newtip becomes
   * the new root.  CAUTION: If newtip is a tip in this case, the resulting
   * tree is degenerate and may not be handled well by other parts of the code.
   * It is therefore recommended that this function be called again immediately
   * with an additional tip node.
   */
{
  node *newfork;
 
  if ( t->root == NULL ) {
    /* TODO: insert single tip */
    return;
  }

  if ( below == NULL ) {
    /* TODO: insert at the root */
    return;
  }
  
  below = t->nodep[below->index - 1];
  newtip = t->nodep[newtip->index-1];

  if ( multf == false ) {
    below = t->nodep[below->index - 1];
    newfork = t->nodep[t->get_fork(t)->index - 1];
    newtip = t->nodep[newtip->index-1];
    if (below->back != NULL)
      below->back->back = newfork;
    newfork->back = below->back;
    below->back = newfork->next->next;
    newfork->next->next->back = below;
    newfork->next->back = newtip;
    newtip->back = newfork->next;
    if (t->root == below)
      t->root = newfork;
  } else {
    newfork = t->get_forknode(t, below->index);
    newfork->next = below->next;
    below->next = newfork;
    hookup(newtip, newfork);
  }
} /* rooted_tree_insert_ */


void buildsimpletree(tree *t, long* enterorder)
{
  node * p = t->nodep[ enterorder[0] - 1];
  node * q = t->nodep[ enterorder[1] - 1];
  node * r = t->nodep[ enterorder[2] - 1];

  hookup(p,q);
  p->v = initialv;
  q->v = initialv;

  t->insert_(t, r, p, false, false);

  t->root = p;

}  /* buildsimpletree */


void rooted_tree_re_move(tree* t, node* item, node** where, boolean doinit)
/* Remove a node from a rooted tree
 *
 * Disconnects item from tree t and if a unifurcation results, joins item's
 * sibling to item's grandparent, freeing item's entire parent fork. If where
 * is given, a pointer to item's former sibling is returned, or NULL if
 * no item could be removed. */
{
  node *whereloc;
  node *p, *q;
  node *fork;
  node *sib;
  if (item == NULL || item->back == NULL) {
    /* TODO Should we die here instead? */
    /* or even set t->root to NULL if item->back == NULL? */
    if (where != NULL)
      *where = NULL;
    return;
  }
    
  if ( count_sibs(item->back) != 2 )  {
    /* removing a node from a multi-furcation is the same in the rooted and unrooted
     * sense */
    generic_tree_re_move(t, item, where, doinit);
    
  } else { /* 2 sibs */
    
    item = t->nodep[item->index-1];
    fork = t->nodep[item->back->index - 1];

    if (item == fork->next->back)
      sib = fork->next->next->back;
    else
      sib = fork->next->back;
    
    if (t->root == fork)
      t->root = sib;
    
    whereloc = sib;
    if ( where ) *where = whereloc;
  
    p = item->back->next->back; /* assumes bifurcation */
    q = item->back->next->next->back;
    if (p != NULL)
      p->back = q;
    if (q != NULL)
      q->back = p;

    t->release_fork(t, fork);
    item->back = NULL;
    if  ( doinit)  {
      inittrav(whereloc);
      inittrav(whereloc->back);
    }
  }
} /* rooted_tree_re_move */


void hsbut(tree* curtree, boolean thorough, boolean jumble, longer seed,
    boolean progress) 
{ /* Heuristic Search for Best Unrooted Tree*/
  long i;
  node* item, *there;
  long *enterorder;
  double bestyet;
  boolean multf;

  enterorder = (long *)Malloc(spp*sizeof(long));
  for (i = 1; i <= spp; i++) 
      enterorder[i - 1] = i;
  if (jumble)
    randumize(seed, enterorder);
  destruct_tree(curtree);	
  buildsimpletree(curtree, enterorder);
  curtree->root = curtree->nodep[enterorder[0] - 1]->back;
  if (progress) {
    printf("Adding species:\n");
    writename(0, 3, enterorder);
    phyFillScreenColor();
    fflush(stdout);
  }
  for (i = 4; i <= spp; i++) {
    bestyet = UNDEFINED;		
    item = curtree->nodep[enterorder[i - 1] - 1];
    there = curtree->root;
    curtree->root = curtree->nodep[enterorder[0] - 1]->back;
    curtree->addtraverse(curtree, item, curtree->root, true, &there, &bestyet,
      NULL, NULL, true, &multf);
    curtree->insert_(curtree, item, there, true, multf);
    curtree->locrearrange(curtree, curtree->nodep[enterorder[0]-1], false, 
        NULL, NULL);
    if (progress) {
      writename(i - 1, 1, enterorder);
      phyFillScreenColor();
      fflush(stdout);
    }
  }
  free(enterorder);
}


void preparetree(tree* t) 
{ /* throw all the forknodes onto the stack so treeread can use them */
  node* p;
  long i;

  while( !Slist_isempty(t->free_forks) ) {
    p = t->get_fork(t);
    t->release_forknode(t, p->next->next);
    t->release_forknode(t, p->next);
    t->release_forknode(t, p);
  }
  for ( i = spp ; i < t->nonodes ; i++ )
    t->nodep[i] = NULL;
}


void fixtree(tree* t)
{ /* after a treeread */
  long i;
  for ( i = spp ; i < t->nonodes ; i++ ) {
    if ( t->nodep[i] == NULL ) {
      t->nodep[i] = t->get_forknode(t, i+1);
      t->nodep[i]->next = t->get_forknode(t, i+1);
      t->nodep[i]->next->next = t->get_forknode(t, i+1);
      t->nodep[i]->next->next->next = t->nodep[i];
      t->release_fork(t, t->nodep[i]);
    }
    else if ( t->nodep[i]->back == NULL && t->nodep[i]->index != t->root->index )
      t->release_fork(t, t->nodep[i]);
  }
}


void arbitrary_resolve(tree* t) 
{ /* gets rid of all multifurcations arbitrarily */
  node *where, *item; 
  long i;

  for ( i = spp ; i < t->nonodes ; i++ ) {
    if ( count_sibs(t->nodep[i]) > 2 ) {
      item = t->nodep[i]->back;
      t->re_move(t, item, &where, false);
      t->insert_(t, item, where, false, false);
      i--; /* do it again, just in case it still multifurcs */
    }
  }
}


void writename(long start, long n, long *enterorder)
{ /* write species name and number in entry order */
  long i, j;

  for (i = start; i < start+n; i++) {
    printf(" %3ld. ", i+1);
    for (j = 0; j < nmlngth; j++)
      putchar(nayme[enterorder[i] - 1][j]);
    putchar('\n');
    fflush(stdout);
  }
}  /* writename */
