from django.db import models


# Create your models here.
class AminoAcid(models.Model):
    """data_storage for the form in the db"""
    X_Progress_ID = models.CharField(max_length=33, blank=True, default="")
    txt = models.FileField(upload_to='AminoAcid/fasta/', null=True, blank=False)
    host = models.CharField(max_length=20, default="413997", blank=True)
    host_threshold = models.FloatField(blank=True, default=0.1)
    max_generation = models.IntegerField(blank=True, default=10)
    # local_homopolymer_threshold = models.IntegerField(blank=True, default=4)
    restriction_enzymes = models.CharField(max_length=200, default="NdeI XhoI HpaI PstI EcoRV NcoI BamHI".split(),
                                           blank=True)  # needs to be split
    splice_sites = models.BooleanField(default=False, blank=True)
    start_sites = models.BooleanField(default=False, blank=True)

    gc_richness_max = models.CharField(max_length=20, blank=True)
    gc_richness_chunk_size = models.CharField(max_length=20, blank=True)
    local_host_profile = models.CharField(max_length=20, blank=True, default=None)


print("model")
