from django import forms
from .models import AminoAcid


class AminoAcidForm(forms.ModelForm):
    class Meta:
        model = AminoAcid
        fields = ("txt",
                  'host',
                  'host_threshold',
                  'max_generation',
                  'restriction_enzymes',
                  'splice_sites',
                  'start_sites',
                  'gc_richness_max',
                  'gc_richness_chunk_size',
                  'local_host_profile',)


print("form")
