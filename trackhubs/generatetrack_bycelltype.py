import glob, os
import trackhub

# First we initialize the components of a track hub

hub, genomes_file, genome, trackdb = trackhub.default_hub(
    hub_name="Bcell_ATACseq",
    short_label='Bcell_ATACseq',
    long_label='Ki67KO in B cells ATACseq',
    genome="mm10",
    email="yan.a@wehi.edu.au")

trackA = trackhub.CompositeTrack(
    shortLabel='PreProB',
    longLabel='Pre Pro B cells',
    visibility='full',
    tracktype='bigWig',
    maxHeightPixels='100:30:10',
    name='A',
)

trackB = trackhub.CompositeTrack(
    shortLabel='ProB',
    longLabel='Pro B cells',
    visibility='full',
    tracktype='bigWig',
    maxHeightPixels='100:30:10',
    name='B',
)

trackC = trackhub.CompositeTrack(
    shortLabel='SmPreB',
    longLabel='Small Pre B cells',
    visibility='full',
    tracktype='bigWig',
    maxHeightPixels='100:30:10',
    name='C',
)

trackdb.add_tracks(trackA)
trackdb.add_tracks(trackB)
trackdb.add_tracks(trackC)

def color_from_filename(fn):
    """
    Figure out a nice color for a track, depending on its filename.
    """
    key = os.path.basename(fn).split('_')[0].lstrip('0123456789')
    colors = {
        'A': '#F8766D',
        'B': '#00BA38',
        'C': '#619CFF',
    } ## define dictionary here
    return trackhub.helpers.hex2rgb(colors[key])

for bigwig in glob.glob('bw/*A_*.bw'):

    name = trackhub.helpers.sanitize(os.path.basename(bigwig))

    track = trackhub.Track(
        name=name,          # track names can't have any spaces or special chars.
        source=bigwig,      # filename to build this track from
        visibility='full',  # shows the full signal
        color=color_from_filename(os.path.basename(bigwig)),  # color_from_filename(os.path.basename(bigwig)), 
        autoScale='group',     # allow the track to autoscale
        smoothingWindow='7',
        tracktype='bigWig', # required when making a track
    )
    trackA.add_subtrack(track)

for bigwig in glob.glob('bw/*B_*.bw'):

    name = trackhub.helpers.sanitize(os.path.basename(bigwig))

    track = trackhub.Track(
        name=name,          # track names can't have any spaces or special chars.
        source=bigwig,      # filename to build this track from
        visibility='full',  # shows the full signal
        color=color_from_filename(os.path.basename(bigwig)),  # 
        autoScale='group',     # allow the track to autoscale
        smoothingWindow='7',
        tracktype='bigWig', # required when making a track
    )
    trackB.add_subtrack(track)

for bigwig in glob.glob('bw/*C_*.bw'):

    name = trackhub.helpers.sanitize(os.path.basename(bigwig))

    track = trackhub.Track(
        name=name,          # track names can't have any spaces or special chars.
        source=bigwig,      # filename to build this track from
        visibility='full',  # shows the full signal
        color=color_from_filename(os.path.basename(bigwig)),  # 
        autoScale='group',     # allow the track to autoscale
        smoothingWindow='7',
        tracktype='bigWig', # required when making a track
    )
    trackC.add_subtrack(track)

trackhub.upload.upload_hub(hub=hub, host='localhost', remote_dir='trackhub')

