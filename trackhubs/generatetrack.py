# import glob, os
# import trackhub

# # First we initialize the components of a track hub

# hub, genomes_file, genome, trackdb = trackhub.default_hub(
#     hub_name="Bcell_ATACseq",
#     short_label='Bcell_ATACseq',
#     long_label='Ki67KO in B cells ATACseq',
#     genome="mm10",
#     email="yan.a@wehi.edu.au")

# trackAWT = trackhub.AggregateTrack(
#     shortLabel='PreProB_WT',
#     longLabel='Pre Pro B cells WT',
#     aggregate='transparentOverlay',
#     visibility='full',
#     tracktype='bigWig',
#     viewLimits='0:250',
#     maxHeightPixels='100:30:10',
#     showSubtrackColorOnUi='on',
#     name='A_WT',
# )

# trackAKO = trackhub.AggregateTrack(
#     shortLabel='PreProB_KO',
#     longLabel='Pre Pro B cells KO',
#     aggregate='transparentOverlay',
#     visibility='full',
#     tracktype='bigWig',
#     viewLimits='0:250',    
#     maxHeightPixels='100:30:10',
#     showSubtrackColorOnUi='on',
#     name='A_KO',
# )

# trackBWT = trackhub.AggregateTrack(
#     shortLabel='ProB_WT',
#     longLabel='Pro B cells WT',
#     aggregate='transparentOverlay',
#     visibility='full',
#     tracktype='bigWig',
#     viewLimits='0:250',
#     maxHeightPixels='100:30:10',
#     showSubtrackColorOnUi='on',
#     name='B_WT',
# )

# trackBKO = trackhub.AggregateTrack(
#     shortLabel='ProB_KO',
#     longLabel='Pro B cells KO',
#     aggregate='transparentOverlay',
#     visibility='full',
#     tracktype='bigWig',
#     viewLimits='0:250',
#     maxHeightPixels='100:30:10',
#     showSubtrackColorOnUi='on',
#     name='B_KO',
# )

# trackCWT = trackhub.AggregateTrack(
#     shortLabel='SmPreB_WT',
#     longLabel='Small Pre B cells WT',
#     aggregate='transparentOverlay',
#     visibility='full',
#     tracktype='bigWig',
#     viewLimits='0:250',
#     maxHeightPixels='100:30:10',
#     showSubtrackColorOnUi='on',
#     name='C_WT',
# )

# trackCKO = trackhub.AggregateTrack(
#     shortLabel='SmPreB_KO',
#     longLabel='Small Pre B cells KO',
#     aggregate='transparentOverlay',
#     visibility='full',
#     tracktype='bigWig',
#     viewLimits='0:250',
#     maxHeightPixels='100:30:10',
#     showSubtrackColorOnUi='on',
#     name='C_KO',
# )

# trackdb.add_tracks(trackAWT)
# trackdb.add_tracks(trackAKO)
# trackdb.add_tracks(trackBWT)
# trackdb.add_tracks(trackBKO)
# trackdb.add_tracks(trackCWT)
# trackdb.add_tracks(trackCKO)

# for bigwig in glob.glob('bw/[1-5]A*.bw'):

#     name = trackhub.helpers.sanitize(os.path.basename(bigwig))

#     track = trackhub.Track(
#         name=name,          # track names can't have any spaces or special chars.
#         source=bigwig,      # filename to build this track from
#         visibility='full',  # shows the full signal
#         color='243,108,101',  # 
#         autoScale='on',     # allow the track to autoscale
#         smoothingWindow='7',
#         tracktype='bigWig', # required when making a track
#     )
#     trackAWT.add_subtrack(track)

# for bigwig in glob.glob('bw/[7-9]A*.bw') + (glob.glob('bw/1[0-1]A*.bw')):

#     name = trackhub.helpers.sanitize(os.path.basename(bigwig))

#     track = trackhub.Track(
#         name=name,          # track names can't have any spaces or special chars.
#         source=bigwig,      # filename to build this track from
#         visibility='full',  # shows the full signal
#         color='243,108,101',  # 
#         autoScale='on',     # allow the track to autoscale
#         smoothingWindow='7',
#         tracktype='bigWig', # required when making a track
#     )
#     trackAKO.add_subtrack(track)

# for bigwig in glob.glob('bw/[1-5]B*.bw'):

#     name = trackhub.helpers.sanitize(os.path.basename(bigwig))

#     track = trackhub.Track(
#         name=name,          # track names can't have any spaces or special chars.
#         source=bigwig,      # filename to build this track from
#         visibility='full',  # shows the full signal
#         color='16,169,71',  # 
#         autoScale='on',     # allow the track to autoscale
#         smoothingWindow='7',
#         tracktype='bigWig', # required when making a track
#     )
#     trackBWT.add_subtrack(track)

# for bigwig in glob.glob('bw/[7-9]B*.bw') + (glob.glob('bw/1[0-1]B*.bw')):

#     name = trackhub.helpers.sanitize(os.path.basename(bigwig))

#     track = trackhub.Track(
#         name=name,          # track names can't have any spaces or special chars.
#         source=bigwig,      # filename to build this track from
#         visibility='full',  # shows the full signal
#         color='16,169,71',  # 
#         autoScale='on',     # allow the track to autoscale
#         smoothingWindow='7',
#         tracktype='bigWig', # required when making a track
#     )
#     trackBKO.add_subtrack(track)

# for bigwig in glob.glob('bw/[1-5]C*.bw'):

#     name = trackhub.helpers.sanitize(os.path.basename(bigwig))

#     track = trackhub.Track(
#         name=name,          # track names can't have any spaces or special chars.
#         source=bigwig,      # filename to build this track from
#         visibility='full',  # shows the full signal
#         color='100,137,196',  # 
#         autoScale='on',     # allow the track to autoscale
#         smoothingWindow='7',
#         tracktype='bigWig', # required when making a track
#     )
#     trackCWT.add_subtrack(track)

# for bigwig in glob.glob('bw/[7-9]C*.bw') + (glob.glob('bw/1[0-1]C*.bw')):

#     name = trackhub.helpers.sanitize(os.path.basename(bigwig))

#     track = trackhub.Track(
#         name=name,          # track names can't have any spaces or special chars.
#         source=bigwig,      # filename to build this track from
#         visibility='full',  # shows the full signal
#         color='100,137,196',  # 
#         autoScale='on',     # allow the track to autoscale
#         smoothingWindow='7',
#         tracktype='bigWig', # required when making a track
#     )
#     trackCKO.add_subtrack(track)

# trackhub.upload.upload_hub(hub=hub, host='localhost', remote_dir='trackhub')

