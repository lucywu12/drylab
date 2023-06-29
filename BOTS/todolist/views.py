from django.shortcuts import render, redirect
from django.core.files.storage import FileSystemStorage
# Create your views here.
from django.shortcuts import render, redirect
from .models import AminoAcid
from .forms import AminoAcidForm
from django.shortcuts import render
from Bio import SeqIO
from django.core.cache import cache
from django.http import HttpResponse, HttpResponseServerError
from todolist.SPEA2_Test.Test_Main import test_main
from todolist.test_functions import is_fasta


def index(request):  # the index view
    print('index')

    def prepareFormForOptimization(form2):
        """Cleans data model into a structure that can easily be transferred through JSON"""
        print("prepareFormForOptimization Start")
        my_set = form2.cleaned_data
        input_file = my_set['txt'].open()
        my_set['txt'] = SeqIO.parse(input_file, "fasta")
        return my_set

    if request.method == 'POST':
        print('POST form')
        form = AminoAcidForm(request.POST, request.FILES)
        if form.is_valid():
            # generate the optimized sequence
            form = prepareFormForOptimization(form)
            print(form)
            if is_fasta(form['txt'].filename):
                test_main()
                return render(request, "index_wrap.html", {'form': form})
            return render(request, "index_wrap.html", {'form': form})

    else:  # not a post, generate empty form
        print('not POST form')
        form = AminoAcidForm(
            initial={"txt": "input.fasta",
                     'host': '413997',
                     'host_threshold': '0.1',
                     'max_generation': '10',
                     'max_relax': '0.1',
                     'restriction_enzymes': 'Ndel Xhol Hpal Pstl EcoRV Ncol BamHI',
                     'splice_sites': False,
                     'start_sites': False,
                     'gc_richness_max': '0.58',
                     'gc_richness_chunk_size': '118',
                     'local_host_profile': None,
                     })
    return render(request, "index_wrap.html", {'form': form})


def index_wrap(request):
    print("index_wrap - view/template")
    return render(request, "index_wrap.html")


def helpFAQ(request):
    print("help - view/template")
    return render(request, "helpFAQ.html")


# A view to report back on upload progress:
def upload_progress(request, form):
    """
    Return JSON object with information about the progress of an upload.
    """
    print(form)
    progress_id = ''
    if 'X_Progress_ID' in request.GET:
        progress_id = request.GET['X_Progress_ID']
    elif 'X_Progress_ID' in request.META:
        progress_id = request.META['X_Progress_ID']
    if progress_id:
        from django.utils import simplejson
        cache_key = "%s_%s" % (request.META['REMOTE_ADDR'], progress_id)
        data = cache.get(cache_key)
        return HttpResponse(simplejson.dumps(data))
    else:
        return HttpResponseServerError('Server Error: You must provide X_Progress_ID header or query param.')


print("views")
