
import copy

# Declare two classes
class Fragment:
    pass

class Fragments:
    pass

class DATAFUNC:
    @staticmethod
    def List(data1, data2):
        if data1 is None or data2 is None:
            return None
        data = []
        if type(data1)==list:
            data += data1
        else:
            data.append(data1)
        if type(data2)==list:
            data += data2
        else:
            data.append(data2)
        return data
    def Max(data1, data2):
        if data1 is None or data2 is None:
            return None
        return max(data1, data2)
    def Min(data1, data2):
        if data1 is None or data2 is None:
            return None
        return min(data1, data2)
    def Mean(data1, data2):
        if data1 is None or data2 is None:
            return None
        return (data1+data2)/2

# Define two classes
# The functions with name end with "_" symbol is mutable functions
# Such as append_, append2_, uniq_, store_as_list_, store_as_dict_,
# TODO: 
# 1. Define copy parameters
# 2. About other parameters of Fragment object

class Fragments:
    """
    This class is used to organize a lot of genomic fragments (regions)
    """
    @staticmethod
    def uniq_fragment_list_(fragment_list: list, data_func=None) -> None:
        """
        Input a sorted fragment list, this function will uniq the fragments, the overlaped regions will be combined
        fragment_list           -- A list of Fragment
        data_func               -- How to process the data field if two Fragment overlaped

        Warning: fragment_list_1 and fragment_list_2 should be sorted

        Return None
        """
        i = len(fragment_list)-1
        while i>0:
            if fragment_list[i]==fragment_list[i-1]:
                # If overlap, then merge
                fragment_list[i-1] = fragment_list[i-1]+fragment_list[i]
                if data_func:
                    # If a data function is provided
                    fragment_list[i-1].data = data_func(fragment_list[i-1].data, fragment_list[i].data)
                del fragment_list[i]
            i -= 1
    @staticmethod
    def split_fragment_list(fragment_list_1: list, fragment_list_2: list, overlap='left') -> Fragments:
        """
        Input two sorted fragment lists, this function will subtract the return the Fragment from fragment_list_1 which not overlap with fragment_list_2
        
        fragment_list_1             -- A list of Fragment
        fragment_list_2             -- A list of Fragment
        overlap                     -- Leave the 'left' or 'right' or 'both' of overlaped regions into both set
        
        Warning: fragment_list_1 and fragment_list_2 should be sorted

        Return Fragments(Spec for 1), Fragments(Overlap), Fragments(Spec for 2)
        """
        assert overlap in ('left', 'both', 'right')
        indicator_1 = [0] * len(fragment_list_1)
        indicator_2 = [0] * len(fragment_list_2)
        j_start = 0
        for i in range(len(indicator_1)):
            f1 = fragment_list_1[i]
            while j_start<len(fragment_list_2) and fragment_list_2[j_start]<f1:
                j_start += 1
            for j in range(j_start, len(fragment_list_2)):
                f2 = fragment_list_2[j]
                if f1==f2:
                    indicator_1[i] = 1
                    indicator_2[j] = 1
                elif f1<f2:
                    break
        frag_1 = Fragments(store='list')   # Spec for 1
        frag_ol = Fragments(store='list')  # Overlap
        frag_2 = Fragments(store='list')   # Spec for 2
        for i,ind in enumerate(indicator_1):
            if ind == 1:
                if overlap in ('left','both'):
                    frag_ol.append_(fragment_list_1[i])
            else:
                frag_1.append_(fragment_list_1[i])
        for j,ind in enumerate(indicator_2):
            if ind == 1:
                if overlap in ('right','both'):
                    frag_ol.append_(fragment_list_2[j])
            else:
                frag_2.append_(fragment_list_2[j])
        return frag_1, frag_ol, frag_2
    def __init__(self, store='list'):
        """
        store           -- Store as list or dict
        """
        assert store in ('list', 'dict')
        
        self.store = store
        
        self.frag_list = []
        self.frag_dict = {}
        
        self.is_sorted = False
    def copy(self):
        import copy
        return copy.deepcopy(self)
    def append_(self, frag: Fragment) -> None:
        """
        Append a object of Fragment
        """
        if self.store == 'list':
            self.frag_list.append( frag )
        elif self.store == 'dict':
            try:
                self.frag_dict[ frag.seqID ].append( frag )
            except KeyError:
                self.frag_dict[ frag.seqID ] = [ frag ]
        self.is_sorted = False
    def append2_(self, seqID, start, end, strand='+', data=None) -> None:
        """
        Append a object of Fragment
        seqID           -- Sequence ID
        start           -- Start coordinate
        end             -- Stop coordinate
        strand          -- + / -
        """
        self.append_( Fragment(seqID, start, end, strand, data) )
    def sort_(self) -> None:
        """
        Sort fragments with coordinate
        """
        if self.is_sorted: return 
        if self.store == 'list':
            self.frag_list.sort()
        elif self.store == 'dict':
            for seqID in self.frag_dict:
                self.frag_dict[seqID].sort()
        self.is_sorted = True
    def uniq_(self, data_func=DATAFUNC.List) -> None:
        """
        Combine overlapped fragments
        """
        self.sort_()
        if self.store == 'list':
            Fragments.uniq_fragment_list_(self.frag_list, data_func=data_func)
        elif self.store == 'dict':
            for seqID in self.frag_dict:
                Fragments.uniq_fragment_list_(self.frag_dict[seqID], data_func=data_func)
    def store_as_list_(self) -> None:
        """
        Store the object as list
        """
        if self.store == 'dict':
            self.frag_list.clear()
            for seqID in self.frag_dict:
                self.frag_list += self.frag_dict[seqID]
            self.store = 'list'
            self.frag_dict.clear()
            self.is_sorted = False
    def store_as_dict_(self) -> None:
        """
        Store the object as dict
        """
        if self.store == 'list':
            self.frag_dict.clear()
            for frag in self.frag_list:
                try:
                    self.frag_dict[frag.seqID].append(frag)
                except KeyError:
                    self.frag_dict[frag.seqID] = [ frag ]
            self.store = 'dict'
            self.frag_list.clear()
    def len(self) -> int:
        """
        Return the number of elements
        """
        if self.store == 'list':
            return len(self.frag_list)
        elif self.store == 'dict':
            return sum( [ len(self.frag_dict[seqID]) for seqID in self.frag_dict ] )
    def __len__(self):
        return self.len()
    def split_fragment(self, other: Fragments, overlap='left') -> Fragments:
        """
        Split two Fragments object into uniq_Fragments_1, Overlaped, uniq_Fragments_2
        
        other           -- A object of Fragments
        overlap         -- Leave the 'left' or 'right' or 'both' of overlaped regions into both set
        
        Return Fragments(Spec for 1), Fragments(Overlap), Fragments(Spec for 2)
        """
        assert overlap in ('left', 'both', 'right')
        #assert self.store == other.store, "Different store type"
        
        self.store_as_list_()
        other.store_as_list_()
        
        self.sort_()
        other.sort_()
        
        return Fragments.split_fragment_list(self.frag_list, other.frag_list, overlap=overlap)
    def remove_(self, other: Fragments) -> None:
        """
        Remove those fragments which overlaped with other.
        """
        frag_1, frag_ol, frag_2 = self.split_fragment(other, overlap='left')
        self.frag_list.clear()
        self.is_sorted = False
        self.frag_list = frag_1.frag_list
    def __add__(self, other: Fragments) -> Fragments:
        """
        Combine two fragments
        other           -- A object of Fragments
        Return a object of Fragments
        """
        #assert self.store == other.store, "Different store type"
        frags = Fragments( store='list' )
        
        self.store_as_list_()
        other.store_as_list_()
        
        frags.frag_list = self.frag_list[:] + other.frag_list[:]
        
        return frags
    def __str__(self):
        return f"<{self.len()} framents store as {self.store}>"
    def __repr__(self):
        return f"<{self.len()} framents store as {self.store}>"
    def show(self, n=10):
        """
        Print the first n elements
        """
        if self.store=='list':
            print( self.frag_list[:n] )
        elif self.store=='dict':
            i = 0
            for seqID in self.frag_dict:
                if i>n: break
                i += 1
                print( seqID )
                print( self.frag_dict[seqID][:n] )

class Fragment:
    def __init__(self, seqID: str, start: int, end: int, strand='+', data=None):
        """
        seqID           -- Sequence ID
        start           -- The start coordinate, closed, 1-based
        end             -- The stop coodinate, open
        strand          -- + / -
        data            -- The data field contains other information
        """
        assert isinstance(start, int)
        assert isinstance(end, int)
        assert start < end
        self.seqID = seqID
        self.start = start
        self.end = end
        self.strand = strand
        self.data = data # the data field 
    def len(self):
        return self.end - self.start
    def copy(self):
        import copy
        return copy.deepcopy(self)
    def equal(self, other):
        """
        The coordinate is identical
        """
        return self.seqID==other.seqID and self.strand==other.strand and self.start==other.start and self.end==other.end
    def __len__(self):
        return self.len()
    def __eq__(self, other):
        """
        If two region overlap, then it equal
        """
        if self.seqID==other.seqID and self.strand==other.strand:
            if self.start<other.end and other.start<self.end:
                return True
        return False
    def __ne__(self, other):
        return not self.__eq__(other)
    def __lt__(self, other):
        #assert self.seqID == other.seqID
        #assert self.strand == other.strand
        if self.seqID<other.seqID:
            return True
        elif self.seqID>other.seqID:
            return False
        if self.strand<other.strand:
            return True
        elif self.strand>other.strand:
            return False
        if self.end <= other.start:
            return True
        else:
            return False
    def __gt__(self, other):
        return not (self<other or self==other)
    def __le__(self, other):
        return self<other or self==other
    def __ge__(self, other):
        return self>other or self==other
    def __sub__(self, other) -> Fragments:
        """
        Remove another fragment from current fragment
        """
        if not self.__eq__(other):
            raise RuntimeError(f"{self} and {other} does not overlap!")
        else:
            frags = Fragments(store='list')
            if other.start<=self.start:
                if other.end>=self.end:
                    return frags
                else:
                    f = self.copy()
                    f.start = other.end
                    frags.append_( f )
            else:
                f = self.copy()
                if other.end>=self.end:
                    f.end = other.start
                    frags.append_( f )
                else:
                    f.end = other.start
                    frags.append_( f )
                    f = self.copy()
                    f.start = other.end
                    frags.append_( f )
            return frags
    def __add__(self, other):
        """
        Merge two frgments
        """
        if not self.__eq__(other):
            raise RuntimeError(f"{self} and {other} does not overlap!")
        else:
            f = self.copy()
            f.start = min( self.start, other.start )
            f.end = max( self.end, other.end )
            return f
    def __str__(self):
        return f"{self.seqID}({self.strand}):{self.start}-{self.end}"
    def __repr__(self):
        return f"{self.seqID}({self.strand}):{self.start}-{self.end}"





