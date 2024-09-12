import random
from faker import Faker

fake = Faker()


def generate_scientific_name():
    """ Produce names like "J.A. Smith" """
    first_name = fake.first_name()
    middle_name = fake.first_name()
    last_name = fake.last_name()
    return f"{first_name[0]}.{middle_name[0]}. {last_name}"


def generate_authors(num_authors=5):
    return ', '.join([generate_scientific_name() for _ in range(num_authors)])

title_start = [
    "The Impact of ",
    "A comparative Study of ",
    "In-Depth Analysis of ",
    "The effect of ",
    "Exploring ",
    "The role of ",
]

research = ['bioinformatics', 'genomics', 'genetics',
            'microbiology', 'personalized medicine', 'virology',
            "protein folding", 'Computational biology']
hotness = ['blockchain', 'AI', 'Space tourism', 'Quantum computing']

for _ in range(100):
    authors = generate_authors(random.randint(2, 4)) + " et al."
    journal = fake.bs().title()
    title = fake.random_element(title_start)
    title += fake.random_element(research) + " and " + fake.random_element(hotness)
    title += " to " + fake.bs()
    print(f"{title} by {authors} et al. {fake.bs().title()} {fake.year()}")

