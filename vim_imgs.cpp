/*
  vim_imags.c
  This file contains all the iamge feature functions.
*/


#include "vim_imgs.h"

//opencv
#ifdef OPENCV_1X 
	#include "cv.h"
	#include "highgui.h"
	#include "cxcore.h"
#else
	#include "opencv2/core/core.hpp"
	#include "opencv2/highgui/highgui.hpp"
	#include "opencv2/calib3d/calib3d.hpp"
	using namespace cv;
#endif



/******************** Functions prototyped K-D-TREE **********************/
/*
	A function to build a k-d tree database from keypoints in an array.
	features - an array of features
	n - the number of features in features
	return:
	Returns the root of a k-d tree built from features or NULL on error.
*/

struct kd_node* kdtree_build( struct feature* features, int n )
{
	struct kd_node* kd_root;

	if( ! features  ||  n <= 0 )
	{
		return NULL;
	}

	kd_root = kd_node_init( features, n );
	expand_kd_node_subtree( kd_root );

	return kd_root;
}

/*
	Finds an image feature's approximate k nearest neighbors in a k-d tree using
	Best Bin First search.

	kd_root - root of an image feature k-d tree
	feat - image feature for whose neighbors to search
	k - number of neighbors to find
	nbrs - pointer to an array in which to store pointers to neighbors
		   in order of increasing descriptor distance
	max_nn_chks - search is cut off after examining this many tree entries

	return:
	Returns the number of neighbors found and stored in nbrs, or -1 on error.
*/
int kdtree_bbf_knn( struct kd_node* kd_root, struct feature* feat, int k,
					struct feature*** nbrs, int max_nn_chks )
{
	struct kd_node* expl;
	struct min_pq* min_pq;
	struct feature* tree_feat, ** _nbrs;
	struct bbf_data* bbf_data;
	int i, t = 0, n = 0;

	if( ! nbrs  ||  ! feat  ||  ! kd_root )
	{
		return -1;
	}

	_nbrs = (struct feature**)calloc( k, sizeof( struct feature* ) );
	min_pq = minpq_init1();
	minpq_insert( min_pq, kd_root, 0 );
	while( min_pq->n > 0  &&  t < max_nn_chks )
	{
		expl = (struct kd_node*)minpq_extract_min( min_pq );
		if( ! expl )
		{
			goto fail;
		}

		expl = explore_to_leaf( expl, feat, min_pq );
		if( ! expl )
		{
			goto fail;
		}

		for( i = 0; i < expl->n; i++ )
		{
			tree_feat = &expl->features[i];
			bbf_data = (struct bbf_data*)malloc( sizeof( struct bbf_data ) );
			if( ! bbf_data )
			{
				goto fail;
			}
			bbf_data->old_data = tree_feat->feature_data;
			bbf_data->d = descr_dist_sq(feat, tree_feat);
			tree_feat->feature_data = bbf_data;
			n += insert_into_nbr_array( tree_feat, _nbrs, n, k );
		}
		t++;
	}

	minpq_release( &min_pq );
	for( i = 0; i < n; i++ )
	{
		bbf_data = (struct bbf_data*)(_nbrs[i]->feature_data);
		_nbrs[i]->feature_data = bbf_data->old_data;
		free( bbf_data );
	}
	*nbrs = _nbrs;
	return n;

fail:
	minpq_release( &min_pq );
	for( i = 0; i < n; i++ )
	{
		bbf_data =(struct bbf_data*) _nbrs[i]->feature_data;
		_nbrs[i]->feature_data = bbf_data->old_data;
		free( bbf_data );
	}
	free( _nbrs );
	*nbrs = NULL;
	return -1;
}

/*
	De-allocates memory held by a k-d tree
	kd_root - pointer to the root of a k-d tree
*/
void kdtree_release( struct kd_node* kd_root )
{
	if( ! kd_root )
		return;
	kdtree_release( kd_root->kd_left );
	kdtree_release( kd_root->kd_right );
	free( kd_root );
}

/************************ Functions prototyped here **************************/
/*
	Initializes a k-d tree node with a set of features.  The node is not
	expanded, and no ordering is imposed on the features.

	features - an array of image features
	n - number of features

	return:
	Returns an unexpanded kd-tree node.
*/
struct kd_node* kd_node_init( struct feature* features, int n )
{
	struct kd_node* kd_node;

	kd_node = (struct kd_node*)malloc( sizeof( struct kd_node ) );
	memset( kd_node, 0, sizeof( struct kd_node ) );
	kd_node->ki = -1;
	kd_node->features = features;
	kd_node->n = n;

	return kd_node;
}

/*
	Recursively expands a specified k-d tree node into a tree whose leaves
	contain one entry each.
	kd_node - an unexpanded node in a k-d tree
*/
void expand_kd_node_subtree( struct kd_node* kd_node )
{
	/* base case: leaf node */
	if( kd_node->n == 1  ||  kd_node->n == 0 )
	{
		kd_node->leaf = 1;
		return;
	}

	assign_part_key( kd_node );
	partition_features( kd_node );

	if( kd_node->kd_left )
		expand_kd_node_subtree( kd_node->kd_left );
	if( kd_node->kd_right )
		expand_kd_node_subtree( kd_node->kd_right );
}

/*
	Determines the descriptor index at which and the value with which to
	partition a k-d tree node's features.
	kd_node - a k-d tree node
*/
void assign_part_key( struct kd_node* kd_node )
{
	struct feature* features;
	double kv, x, mean, var, var_max = 0;
	double* tmp;
	int d, n, i, j, ki = 0;

	features = kd_node->features;
	n = kd_node->n;
	d = features[0].d;

	/* partition key index is that along which descriptors have most variance */
	for( j = 0; j < d; j++ )
	{
		mean = var = 0;
		for( i = 0; i < n; i++ )
			mean += features[i].descr[j];
		mean /= n;
		for( i = 0; i < n; i++ )
		{
			x = features[i].descr[j] - mean;
			var += x * x;
		}
		var /= n;

		if( var > var_max )
		{
			ki = j;
			var_max = var;
		}
	}

	/* partition key value is median of descriptor values at ki */
	tmp = (double*)calloc( n, sizeof( double ) );
	for( i = 0; i < n; i++ )
		tmp[i] = features[i].descr[ki];
	kv = median_select( tmp, n );
	free( tmp );

	kd_node->ki = ki;
	kd_node->kv = kv;
}

/*
	Finds the median value of an array.  The array's elements are re-ordered
	by this function.

	array - an array; the order of its elelemts is reordered
	n - number of elements in array

	return:
	Returns the median value of array.
*/
double median_select( double* array, int n )
{
	return rank_select( array, n, (n - 1) / 2 );
}

/*
	Finds the element of a specified rank in an array using the linear time
	median-of-medians algorithm by Blum, Floyd, Pratt, Rivest, and Tarjan.
	The elements of the array are re-ordered by this function.

	array - an array; the order of its elelemts is reordered
	n - number of elements in array
	r - the zero-based rank of the element to be selected

	return:
	Returns the element from array with zero-based rank r.
*/
double rank_select( double* array, int n, int r )
{
	double* tmp, med;
	int gr_5, gr_tot, rem_elts, i, j;

	/* base case */
	if( n == 1 )
		return array[0];

	/* divide array into groups of 5 and sort them */
	gr_5 = n / 5;
	gr_tot = cvCeil( n / 5.0 );
	rem_elts = n % 5;
	tmp = array;
	for( i = 0; i < gr_5; i++ )
	{
		insertion_sort( tmp, 5 );
		tmp += 5;
	}
	insertion_sort( tmp, rem_elts );

	/* recursively find the median of the medians of the groups of 5 */
	tmp = (double*)calloc( gr_tot, sizeof( double ) );
	for( i = 0, j = 2; i < gr_5; i++, j += 5 )
		tmp[i] = array[j];
	if( rem_elts )
		tmp[i++] = array[n - 1 - rem_elts/2];
	med = rank_select( tmp, i, ( i - 1 ) / 2 );
	free( tmp );

	/* partition around median of medians and recursively select if necessary */
	j = partition_array( array, n, med );
	if( r == j )
		return med;
	else if( r < j )
		return rank_select( array, j, r );
	else
	{
		array += j+1;
		return rank_select( array, ( n - j - 1 ), ( r - j - 1 ) );
	}
}

/*
	Sorts an array in place into increasing order using insertion sort.

	array - an array
	n - number of elements
*/
void insertion_sort( double* array, int n )
{
	double k;
	int i, j;

	for( i = 1; i < n; i++ )
	{
		k = array[i];
		j = i-1;
		while( j >= 0  &&  array[j] > k )
		{
			array[j+1] = array[j];
			j -= 1;
		}
		array[j+1] = k;
	}
}

/*
	Partitions an array around a specified value.

	array - an array
	n - number of elements
	pivot - value around which to partition

	return:
	Returns index of the pivot after partitioning
*/
int partition_array( double* array, int n, double pivot )
{
	double tmp;
	int p, i, j;

	i = -1;
	for( j = 0; j < n; j++ )
	{
		if( array[j] <= pivot )
		{
			tmp = array[++i];
			array[i] = array[j];
			array[j] = tmp;
			if( array[i] == pivot )
				p = i;
		}
	}
	array[p] = array[i];
	array[i] = pivot;

	return i;
}

/*
	Partitions the features at a specified k-d tree node to create its two
	children.

	kd_node - a k-d tree node whose partition key is set
*/
void partition_features( struct kd_node* kd_node )
{
	struct feature* features, tmp;
	double kv;
	int n, ki, p, i, j = -1;

	features = kd_node->features;
	n = kd_node->n;
	ki = kd_node->ki;
	kv = kd_node->kv;
	for( i = 0; i < n; i++ )
	{
		if( features[i].descr[ki] <= kv )
		{
			tmp = features[++j];
			features[j] = features[i];
			features[i] = tmp;
			if( features[j].descr[ki] == kv )
				p = j;
		}
	}

	tmp = features[p];
	features[p] = features[j];
	features[j] = tmp;

	/* if all records fall on same side of partition, make node a leaf */
	if( j == n - 1 )
	{
		kd_node->leaf = 1;
		return;
	}

	kd_node->kd_left = kd_node_init( features, j + 1 );
	kd_node->kd_right = kd_node_init( features + ( j + 1 ), ( n - j - 1 ) );
}

/*
	Explores a k-d tree from a given node to a leaf.  

	kd_node - root of the subtree to be explored
	feat - feature upon which branching decisions are based
	min_pq - a minimizing priority queue into which tree nodes are placed
		as described above

	return:
	Returns a pointer to the leaf node at which exploration ends or
		NULL on error.
*/
struct kd_node* explore_to_leaf( struct kd_node* kd_node, struct feature* feat,
								struct min_pq* min_pq )
{
	struct kd_node* unexpl, * expl = kd_node;
	double kv;
	int ki;

	while( expl  &&  ! expl->leaf )
	{
		ki = expl->ki;
		kv = expl->kv;

		if( ki >= feat->d )
		{
			return NULL;
		}
		if( feat->descr[ki] <= kv )
		{
			unexpl = expl->kd_right;
			expl = expl->kd_left;
		}
		else
		{
			unexpl = expl->kd_left;
			expl = expl->kd_right;
		}

		if( minpq_insert( min_pq, unexpl, (int) (ABS( kv - feat->descr[ki] ) ) ))
		{
			return NULL;
		}
	}

	return expl;
}

/*
	Inserts a feature into the nearest-neighbor array so that the array remains
	in order of increasing descriptor distance from the search feature.

	feat - feature to be inderted into the array; it's feature_data field
		should be a pointer to a bbf_data with d equal to the squared descriptor
		distance between feat and the search feature
	nbrs- array of nearest neighbors neighbors
	n - number of elements already in nbrs and
	k - maximum number of elements in nbrs

	return:
	If feat was successfully inserted into nbrs, returns 1; otherwise
		returns 0.
*/
int insert_into_nbr_array( struct feature* feat, struct feature** nbrs,
						  int n, int k )
{
	struct bbf_data* fdata, * ndata;
	double dn, df;
	int i, ret = 0;

	if( n == 0 )
	{
		nbrs[0] = feat;
		return 1;
	}

	/* check at end of array */
	fdata = (struct bbf_data*)feat->feature_data;
	df = fdata->d;
	ndata = (struct bbf_data*)nbrs[n-1]->feature_data;
	dn = ndata->d;
	if( df >= dn )
	{
		if( n == k )
		{
			feat->feature_data = fdata->old_data;
			free( fdata );
			return 0;
		}
		nbrs[n] = feat;
		return 1;
	}

	/* find the right place in the array */
	if( n < k )
	{
		nbrs[n] = nbrs[n-1];
		ret = 1;
	}
	else
	{
		nbrs[n-1]->feature_data = ndata->old_data;
		free( ndata );
	}
	i = n-2;
	while( i >= 0 )
	{
		ndata = (struct bbf_data*)nbrs[i]->feature_data;
		dn = ndata->d;
		if( dn <= df )
			break;
		nbrs[i+1] = nbrs[i];
		i--;
	}
	i++;
	nbrs[i] = feat;

	return ret;
}

/*
Creates a new minimizing priority queue.
*/
struct min_pq* minpq_init1()
{
	struct min_pq* min_pq;

	min_pq =(struct min_pq*) malloc( sizeof( struct min_pq ) );
	min_pq->pq_array = (struct pq_node*)calloc( MINPQ_INIT_NALLOCD, sizeof( struct pq_node ) );
	min_pq->nallocd = MINPQ_INIT_NALLOCD;
	min_pq->n = 0;

	return min_pq;
}


/**
	Inserts an element into a minimizing priority queue.

	min_pq - a minimizing priority queue
	data  - the data to be inserted
	key  - the key to be associated with \a data

	return 
	Returns 0 on success or 1 on failure.
*/
int minpq_insert( struct min_pq* min_pq, void* data, int key )
{
	int n = min_pq->n;

	/* double array allocation if necessary */
	if( min_pq->nallocd == n )
	{
		min_pq->nallocd = array_double1( (void**)&min_pq->pq_array, min_pq->nallocd,
										sizeof( struct pq_node ) );
		if( ! min_pq->nallocd )
		{
			return 1;
		}
	}

	min_pq->pq_array[n].data = data;
	min_pq->pq_array[n].key = INT_MAX;
	decrease_pq_node_key( min_pq->pq_array, min_pq->n, key );
	min_pq->n++;

	return 0;
}

/*
	Removes and returns the element of a minimizing priority queue with the
	smallest key.

	min_pq a minimizing priority queue

	return:
	Returns the element of \a min_pq with the smallest key 

*/
void* minpq_extract_min( struct min_pq* min_pq )
{
	void* data;

	if( min_pq->n < 1 )
	{
		return NULL;
	}
	data = min_pq->pq_array[0].data;
	min_pq->n--;
	min_pq->pq_array[0] = min_pq->pq_array[min_pq->n];
	restore_minpq_order( min_pq->pq_array, 0, min_pq->n );

	return data;
}

/*
	De-allocates the memory held by a minimizing priorioty queue
	min_pq - pointer to a minimizing priority queue
*/
void minpq_release( struct min_pq** min_pq )
{
	if( ! min_pq )
	{
		return;
	}
	if( *min_pq  &&  (*min_pq)->pq_array )
	{
		free( (*min_pq)->pq_array );
		free( *min_pq );
		*min_pq = NULL;
	}
}


/************************ Functions prototyped here **************************/
/**
	Decrease a minimizing pq element's key, rearranging the pq if necessary

	pq_array -  minimizing priority queue array
	i -   index of the element whose key is to be decreased
	key - new value of element <EM>i</EM>'s key; if greater than current
		key, no action is taken
*/
void decrease_pq_node_key( struct pq_node* pq_array, int i, int key )
{
	struct pq_node tmp;

	if( key > pq_array[i].key )
		return;

	pq_array[i].key = key;
	while( i > 0  &&  pq_array[i].key < pq_array[parent(i)].key )
	{
		tmp = pq_array[parent(i)];
		pq_array[parent(i)] = pq_array[i];
		pq_array[i] = tmp;
		i = parent(i);
	}
}

/*
	Recursively restores correct priority queue order to a minimizing pq array

	pq_array - a minimizing priority queue array
	i -  index at which to start reordering
	n -  number of elements in \a pq_array
*/
void restore_minpq_order( struct pq_node* pq_array, int i, int n )
{
	struct pq_node tmp;
	int l, r, min = i;

	l = left( i );
	r = right( i );
	if( l < n )
		if( pq_array[l].key < pq_array[i].key )
			min = l;
	if( r < n )
		if( pq_array[r].key < pq_array[min].key )
			min = r;

	if( min != i )
	{
		tmp = pq_array[min];
		pq_array[min] = pq_array[i];
		pq_array[i] = tmp;
		restore_minpq_order( pq_array, min, n );
	}
}

/*************************** Function Definitions ****************************/

void fatal_error(char* format, ...)
{
	va_list ap;
	fprintf( stderr, "Error: ");

	va_start( ap, format );
	vfprintf( stderr, format, ap );
	
	va_end( ap );
	fprintf( stderr, "\n" );
	abort();
}


/*
	Doubles the size of an array with error checking 
	array - pointer to an array whose size is to be doubled
	n	  - number of elements allocated for \a array
	size  - size in bytes of elements in \a array

	return : 
	Returns the new number of elements allocated for \a array.  
	If no memory is available, returns 0 and frees array.
*/
int array_double1( void** array, int n, int size )
{
	void* tmp;

	tmp = realloc( *array, 2 * n * size );
	if( ! tmp )
	{
		if( *array )
			free( *array );
		*array = NULL;
		return 0;
	}
	*array = tmp;
	return n*2;
}


/************************************************************************/
/************************************************************************/
/************************************************************************/

/*
	Allocates and initializes a new feature
	return:
	 Returns a pointer to the new feature
*/
struct feature* new_feature( void )
{
	struct feature* feat;
	struct detection_data* ddata;

	feat = (struct feature*)malloc( sizeof( struct feature ) );
	memset( feat, 0, sizeof( struct feature ) );
	ddata = (struct detection_data*)malloc( sizeof( struct detection_data ) );
	memset( ddata, 0, sizeof( struct detection_data ) );
	feat->feature_data = ddata;

	return feat;
}

/*
	Gaussian smooths an orientation histogram.

	hist - an orientation histogram
	n - number of bins
*/
void smooth_ori_hist( double* hist, int n )
{
	double prev, tmp, h0 = hist[0];
	int i;

	prev = hist[n-1];
	for( i = 0; i < n; i++ )
	{
		tmp = hist[i];
		hist[i] = 0.25 * prev + 0.5 * hist[i] + 
			0.25 * ( ( i+1 == n )? h0 : hist[i+1] );
		prev = tmp;
	}
}

/*
	Finds the magnitude of the dominant orientation in a histogram

	hist - an orientation histogram
	n - number of bins

	return 
	Returns the value of the largest bin in hist
*/
double dominant_ori( double* hist, int n )
{
	double omax;
	int maxbin, i;

	omax = hist[0];
	maxbin = 0;
	for( i = 1; i < n; i++ )
		if( hist[i] > omax )
		{
			omax = hist[i];
			maxbin = i;
		}
	return omax;
}


/*
	Makes a deep copy of a feature
	feat - feature to be cloned
	 return:
	 Returns a deep copy of feat
*/
struct feature* clone_feature( struct feature* feat )
{
	struct feature* new_feat;
	struct detection_data* ddata;

	new_feat = new_feature();
	ddata = feat_detection_data( new_feat );
	memcpy( new_feat, feat, sizeof( struct feature ) );
	memcpy( ddata, feat_detection_data(feat), sizeof( struct detection_data ) );
	new_feat->feature_data = ddata;

	return new_feat;
}

/*
	Interpolates an entry into the array of orientation histograms that form
	the feature descriptor.

	hist  - 2D array of orientation histograms
	rbin  - sub-bin row coordinate of entry
	cbin  - sub-bin column coordinate of entry
	obin  - sub-bin orientation coordinate of entry
	mag  - size of entry
	d  - width of 2D array of orientation histograms
	n  - number of bins per orientation histogram
*/
void interp_hist_entry( double*** hist, double rbin, double cbin,
					   double obin, double mag, int d, int n )
{
	double d_r, d_c, d_o, v_r, v_c, v_o;
	double** row, * h;
	int r0, c0, o0, rb, cb, ob, r, c, o;

	r0 = cvFloor( rbin );
	c0 = cvFloor( cbin );
	o0 = cvFloor( obin );
	d_r = rbin - r0;
	d_c = cbin - c0;
	d_o = obin - o0;

	/*
	The entry is distributed into up to 8 bins.  Each entry into a bin
	is multiplied by a weight of 1 - d for each dimension, where d is the
	distance from the center value of the bin measured in bin units.
	*/
	for( r = 0; r <= 1; r++ )
	{
		rb = r0 + r;
		if( rb >= 0  &&  rb < d )
		{
			v_r = mag * ( ( r == 0 )? 1.0 - d_r : d_r );
			row = hist[rb];
			for( c = 0; c <= 1; c++ )
			{
				cb = c0 + c;
				if( cb >= 0  &&  cb < d )
				{
					v_c = v_r * ( ( c == 0 )? 1.0 - d_c : d_c );
					h = row[cb];
					for( o = 0; o <= 1; o++ )
					{
						ob = ( o0 + o ) % n;
						v_o = v_c * ( ( o == 0 )? 1.0 - d_o : d_o );
						h[ob] += v_o;
					}
				}
			}
		}
	}
}

/*
	Converts the 2D array of orientation histograms into a feature's descriptor
	vector.

	hist  - 2D array of orientation histograms
	d  - width of hist
	n  - bins per histogram
	feat - feature into which to store descriptor
*/
void hist_to_descr_descriptor_vector( double*** hist, int d, int n, struct feature* feat )
{
	int int_val, i, r, c, o, k = 0;

	for( r = 0; r < d; r++ )
		for( c = 0; c < d; c++ )
			for( o = 0; o < n; o++ )
				feat->descr[k++] = hist[r][c][o];

	feat->d = k;
	normalize_descriptor_vector( feat );
	for( i = 0; i < k; i++ )
		if( feat->descr[i] > SIFT_DESCR_MAG_THR )
			feat->descr[i] = SIFT_DESCR_MAG_THR;
	normalize_descriptor_vector( feat );

	/* convert floating-point descriptor to integer valued descriptor */
	for( i = 0; i < k; i++ )
	{
		int_val = (int) (SIFT_INT_DESCR_FCTR * feat->descr[i]);
		feat->descr[i] = MIN( 255, int_val );
	}
}

/*
	Normalizes a feature's descriptor vector to unitl length
	feat - feature
*/
void normalize_descriptor_vector( struct feature* feat )
{
	double cur, len_inv, len_sq = 0.0;
	int i, d = feat->d;

	for( i = 0; i < d; i++ )
	{
		cur = feat->descr[i];
		len_sq += cur*cur;
	}
	len_inv = 1.0 / sqrt( len_sq );
	for( i = 0; i < d; i++ )
		feat->descr[i] *= len_inv;
}

/*
	Compares features for a decreasing-scale ordering.  
	Intended for use with CvSeqSort

	feat1 -  first feature
	feat2 -  second feature

	return:
	Returns 1 if feat1's scale is greater than feat2's, 
		   -1 if vice versa,
			0 if their scales are equal
*/
int feature_cmpare( void* feat1, void* feat2)
{
	struct feature* f1 = (struct feature*) feat1;
	struct feature* f2 = (struct feature*) feat2;

	if( f1->scl < f2->scl )
		return 1;
	if( f1->scl > f2->scl )
		return -1;
	return 0;
}

/*
	De-allocates memory held by a descriptor histogram

	hist -  pointer to a 2D array of orientation histograms
	d - width of hist
*/
void release_descr_hist( double**** hist, int d )
{
	int i, j;

	for( i = 0; i < d; i++)
	{
		for( j = 0; j < d; j++ )
			free( (*hist)[i][j] );
		free( (*hist)[i] );
	}
	free( *hist );
	*hist = NULL;
}

//计算矢量的欧式距离
//////////////////////////////////////////////////////////////////////////
double descr_dist_sq( struct feature* f1, struct feature* f2 )
{
	double diff, dsq = 0;
	float* descr1, * descr2;
	int i, d;	
	
	d = f1->d;
	if( f2->d != d )
		return DBL_MAX;

	descr1 = f1->descr;
	descr2 = f2->descr;
	
	for( i = 0; i < d; i++ )
	{
		diff = descr1[i] - descr2[i];
		dsq += diff*diff;
	}
	return dsq;
}
