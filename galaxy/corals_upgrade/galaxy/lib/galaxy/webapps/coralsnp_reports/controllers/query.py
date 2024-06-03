""" Mixin to help build advanced queries coralsnp for reports interface.
"""
import sqlalchemy as sa


class ReportQueryBuilder:
    def group_by_month(self, column):
        return [sa.func.date_trunc("month", sa.func.date(column))]

    def select_month(self, column):
        return sa.func.date_trunc("month", sa.func.date(column))

    def group_by_day(self, column):
        return [sa.func.date_trunc("day", sa.func.date(column))]

    def select_day(self, column):
        return sa.func.date_trunc("day", sa.func.date(column))
