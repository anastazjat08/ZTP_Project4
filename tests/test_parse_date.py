from literature.pubmed_fetch import parse_date

def test_parse_date_full():
    pub_date = {
        'Year': '2021',
        'Month': 'Jan',
        'Day': '10'
    }

    result = parse_date(pub_date)

    assert result["Year"] == 2021
    assert result["Month"] == "Jan"
    assert result["Day"] == 10

def test_parse_date_partial():
    pub_date = {
        'Year': '2021',
        'Month': 'Jan'
    }

    result = parse_date(pub_date)

    assert result["Year"] == 2021
    assert result["Month"] == "Jan"

def test_parse_date_none():
    result = parse_date({})

    assert result["Year"] is None
    assert result["Month"] is None
    assert result["Day"] is None