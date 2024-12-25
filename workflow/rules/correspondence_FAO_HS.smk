rule correspondence_FAO_HS:
    output:
        'resources/raw_data/correspondence_FAO_HS.json'
    threads: 1
    script:
        '../scripts/correspondence_FAO_HS.py'